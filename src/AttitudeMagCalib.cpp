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
#include <Eigen/Cholesky>
#include <cmath>

using namespace Eigen;

namespace kr {

AttitudeMagCalib::AttitudeMagCalib() { reset(); }

void AttitudeMagCalib::reset() {
  binH_.clear();
  binV_.clear();
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
  
  //  determine local angle
  if (std::abs(localG[2]) < 0.1) {
    //  world vertical is approx. in the local X/Y plane
    const kr::vec3d worldZ = att.matrix() * kr::vec3d(0,0,1);
    const double ang = std::atan2(worldZ[1], worldZ[0]); 
    const int key = (ang + M_PI) / (2*M_PI) * kBinMaxCount;
    binV_[key] = bin;
  } else if (std::abs(localG[2]) > 0.9) {
    //  world vertical is approx. vertical
    const kr::vec3d worldX = att.matrix() * kr::vec3d(1,0,0);
    const double ang = std::atan2(worldX[1], worldX[0]);
    const int key = (ang + M_PI) / (2*M_PI) * kBinMaxCount;
    binH_[key] = bin;
  }
}

bool AttitudeMagCalib::isReady() const {
  if (binV_.size() < kBinMaxCount*8/10) {
    return false;
  }
  if (binH_.size() < kBinMaxCount*8/10) {
    return false;
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
    for (const std::pair<int, SampleBin>& s : binH_) {
      //printf("%f, %f, %f\n", s.second.field[0], s.second.field[1], s.second.field[2]);
      meas.push_back(s.second.field);
    }
    for (const std::pair<int, SampleBin>&s : binV_) {
      //printf("%f, %f, %f\n", s.second.field[0], s.second.field[1], s.second.field[2]);
      meas.push_back(s.second.field);
    }
    const size_t N = meas.size();

    //  fit to sphere
    MatrixXd A(N, 4);
    VectorXd b(N, 1);
    for (size_t i=0; i < meas.size(); i++) {
      A(i,0) = 2*meas[i][0];
      A(i,1) = 2*meas[i][1];
      A(i,2) = 2*meas[i][2];
      A(i,3) = 1;
      b(i,0) = meas[i][0]*meas[i][0] + meas[i][1]*meas[i][1] + 
          meas[i][2]*meas[i][2];
    }
    //  solve system for center
    const Vector4d x = A.colPivHouseholderQr().solve(b);
    bias = x.block<3,1>(0,0);

    double mean_rad = 0.0, mean_rad_sqr = 0.0;
    //  calculate estimate of mean radius
    for (const Vector3d &v : meas) {
      const double r = (v - bias).norm();
      mean_rad += r;
      mean_rad_sqr += r * r;
    }
    mean_rad /= N;
    mean_rad_sqr /= N;

    //  refine with GN-NLS
    Matrix<double, Eigen::Dynamic, 6> J(N, 6);
    Matrix<double, Eigen::Dynamic, 1> r(N, 1);
    MatrixXd W(N,N);
    W.setZero();

    for (int iter = 0; iter < 20; iter++) {
      for (size_t i = 0; i < meas.size(); i++) {
        const double x = (meas[i][0] - bias[0]) / scl[0];
        const double y = (meas[i][1] - bias[1]) / scl[1];
        const double z = (meas[i][2] - bias[2]) / scl[2];
        const double r2 = x*x + y*y + z*z;
        
        //  residual
        r(i,0) = mean_rad*mean_rad - r2;
        //  jacobian
        J(i,0) = -2 * x * x / scl[0];
        J(i,1) = -2 * y * y / scl[1];
        J(i,2) = -2 * z * z / scl[2];
        J(i,3) = -2 * x / scl[0];
        J(i,4) = -2 * y / scl[1];
        J(i,5) = -2 * z / scl[2];
      }

      //  calculate mean squared error
      double sigmaSquared=0;
      for (int j=0; j < r.rows(); j++) {
        sigmaSquared += r(j,0)*r(j,0);
      }
      sigmaSquared /= r.rows();
      //  calculate cauchy weights
      for (int j=0; j < r.rows(); j++) {
        const double errSqr = r(j,0)*r(j,0);
        W(j,j) = 1.0 / (1 + errSqr/sigmaSquared);
      }
      
      Matrix<double, 6, 6> H = J.transpose() * W * J;
      for (int i = 0; i < 3; i++) {  //  prior on scale
        H(i,i) *= 1.01;
      }
      Eigen::LDLT<Matrix<double,6,6>> LDLT(H);
      const Matrix<double, 6, 1> update = LDLT.solve(J.transpose() * W * r);

      scl += update.block<3, 1>(0, 0);
      bias += update.block<3, 1>(3, 0);
    }

    bias_ = bias;
    scale_ = scl;
  } else {
    bias_.setZero();
    scale_.setOnes();
  }

  //  calculate the magnetic reference vector using the yaw samples
  double hAvg = 0.0, vAvg = 0.0;
  
  for (const std::pair<int,SampleBin>& p : binH_) {
    const kr::mat3d wRb = p.second.q.matrix();
    const kr::vec3d rpy = kr::getRPY(wRb);
    const kr::mat3d tilt = kr::rotation_y(rpy[1]) * kr::rotation_x(rpy[0]);
    kr::vec3d level = tilt * p.second.field;
        
    level -= bias_;
    for (int i=0; i < 3; i++) {
      level[i] /= scale_[i];
    }
    
    hAvg += std::sqrt(level[0]*level[0] + level[1]*level[1]);
    vAvg += level[2];
  }
  
  ref_[0] = 0;
  ref_[1] = hAvg / binH_.size(); //  Y = North
  ref_[2] = -std::abs(vAvg) / binH_.size();
  calibrated_ = true;
}

} //  namespace kr
