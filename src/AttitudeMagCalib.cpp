/*
 * AttitudeMagCalib.cpp
 *
 *  Copyright (c) 2013 Gareth Cross. All rights reserved.
 *
 *  This file is part of kr_attitude_eskf.
 *
 *	Created on: 23/06/2014
 *		  Author: gareth
 */

#include <kr_attitude_eskf/AttitudeMagCalib.hpp>
#include <kr_math/SO3.hpp>
#include <cmath>

using namespace Eigen;

namespace kr {

AttitudeMagCalib::AttitudeMagCalib() { reset(); }

void AttitudeMagCalib::reset() {
  for (int i = 0; i < 3; i++) {
    axes_[i].clear();
  }
  calibrated_ = false;
  bias_.setZero();
  scale_.setOnes();
}

void AttitudeMagCalib::appendSample(const kr::quatd &att,
                                    const Vector3d &field) {
  SampleBin bin;
  bin.field = field;
  bin.q = att;
  
  const kr::vec3d localG = att.conjugate().matrix() * kr::vec3d(0,0,1);
  
  if (std::abs(localG[2]) < 0.2) {
    //  world vertical is approx. in the local X/Y plane
    
  } else if (std::abs(localG[2]) > 0.8) {
    //  world vertical is approx. vertical
    
  }
  
  const Vector3d angs = kr::getRPY(att.matrix());

  //  map to bins
  int n_roll = std::floor((angs[0] + M_PI) / (2 * M_PI) * kBinMaxCount);
  int n_pitch = std::floor((angs[1] + M_PI / 2) / M_PI * kBinMaxCount);
  int n_yaw = std::floor((angs[2] + M_PI) / (2 * M_PI) * kBinMaxCount);

  if (axes_[0].find(n_roll) == axes_[0].end() || 
      axes_[1].find(n_pitch) == axes_[1].end() ||
      axes_[2].find(n_yaw) == axes_[2].end()) {
    printf("%f, %f, %f\n", field[0], field[1], field[2]);
  }
  
  axes_[0][n_roll] = bin;
  axes_[1][n_pitch] = bin;
  axes_[2][n_yaw] = bin;
}

bool AttitudeMagCalib::isReady() const {
  //  must collect 90% of samples to be ready
  for (int i = 0; i < 3; i++) {
    if (axes_[i].size() < kBinMaxCount * 9 / 10) {
      return false;
    }
  }
  return true;
}

bool AttitudeMagCalib::isCalibrated() const { return calibrated_; }

void AttitudeMagCalib::calibrate(AttitudeMagCalib::CalibrationType type) {

  if (!isReady()) {
    throw insufficient_data();
  }

  if (type == AttitudeMagCalib::FullCalibration) {
    //  perform full estimation of bias and scale
    Vector3d bias, scl = Vector3d::Ones();

    std::vector<Vector3d> meas;
    for (int i = 0; i < 3; i++) {
      for (const auto &s : axes_[i]) {
        meas.push_back(s.second.field);
      }
    }
    const size_t N = meas.size();

    Vector3d max, min, mean;
    min.setConstant(3, std::numeric_limits<double>::infinity());
    max = -min;
    mean.setZero();

    double mean_rad = 0.0, mean_rad_sqr = 0.0;

    for (const Vector3d &v : meas) {
      const double r = v.norm();
      mean_rad += r;
      mean_rad_sqr += r * r;

      for (int i = 0; i < 3; i++) {
        max[i] = std::max(max[i], v[i]);
        min[i] = std::min(min[i], v[i]);
      }

      mean += v;
    }
    mean /= N;
    mean_rad /= N;
    mean_rad_sqr /= N;

    const float var = mean_rad_sqr - mean_rad * mean_rad;
    const float den = std::sqrt(2 * M_PI * var);

    //  initial estimate of the bias
    bias = (max + min) * 0.25 + mean * 0.5;

    //  refine with GN-NLS
    Matrix<double, Eigen::Dynamic, 6> J(N, 6);
    Matrix<double, Eigen::Dynamic, 1> r(N, 1);
    MatrixXd W;

    //  calculate weights
    W.resize(N, N);
    W.setZero();
    for (size_t i = 0; i < N; i++) {
      const double r = meas[i].norm();
      W(i, i) = std::exp(-0.5 * (r - mean_rad) * (r - mean_rad) / var) / den;
    }

    for (int iter = 0; iter < 10; iter++) {
      for (size_t i = 0; i < meas.size(); i++) {
        const double x = (meas[i][0] - bias[0]) / scl[0];
        const double y = (meas[i][1] - bias[1]) / scl[1];
        const double z = (meas[i][2] - bias[2]) / scl[2];
        const double r2 = x * x + y * y + z * z;

        r[i] = W(i, i) * (mean_rad * mean_rad - r2); //  weighted residual

        //  jacobian
        J(i, 0) = -2 * x * x / scl[0];
        J(i, 1) = -2 * y * y / scl[1];
        J(i, 2) = -2 * z * z / scl[2];
        J(i, 3) = -2 * x / scl[0];
        J(i, 4) = -2 * y / scl[1];
        J(i, 5) = -2 * z / scl[2];
      }

      Matrix<double, 6, 6> H = J.transpose() * W * J;
      for (int i = 0; i < H.rows(); i++) {
        H(i, i) *= 1.01;
      }
      auto LU = H.fullPivLu();

      if (!LU.isInvertible()) {
        throw singular_hessian();
      }

      decltype(H) Hinv = LU.inverse();
      Matrix<double, 6, 1> update = Hinv * J.transpose() * r;

      bias += update.block<3, 1>(3, 0);
      scl += update.block<3, 1>(0, 0);
    }

    bias_ = bias;
    scale_ = scl;
  } else {
    bias_.setZero();
    scale_.setOnes();
  }

  //  calculate the magnetic reference vector using the yaw samples
  double hAvg = 0.0, vAvg = 0.0;

  for (auto i = axes_[2].begin(); i != axes_[2].end(); i++) {
    //  tilt-compensate the components
    const Matrix3d wRb = i->second.q.matrix();
    const Vector3d rpy = kr::getRPY(wRb);
    const Matrix3d tilt = kr::rotation_y(rpy[1]) * kr::rotation_x(rpy[0]);

    Vector3d level = tilt * i->second.field;
    for (int j = 0; j < 3; j++) {
      level[j] = (level[j] - bias_[j]) / scale_[j];
    }

    hAvg += std::sqrt(level[0] * level[0] + level[1] * level[1]);
    vAvg += level[2];
  }

  ref_[0] = 0;
  ref_[1] = hAvg / axes_[2].size(); //  Y = North
  ref_[2] = vAvg / axes_[2].size();
  calibrated_ = true;
}

} //  namespace kr
