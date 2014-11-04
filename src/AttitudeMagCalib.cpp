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

void AttitudeMagCalib::appendSample(const quat &att,
                                    const vec3 &field) {
  SampleBin bin;
  bin.field = field;
  bin.q = att;
  
  const vec3 localG = att.conjugate().matrix() * vec3(0,0,1);
  
  //  determine local angle
  if (std::abs(localG[2]) < 0.1) {
    //  world vertical is approx. in the local X/Y plane
    const vec3 worldZ = att.matrix() * vec3(0,0,1);
    const scalar_t ang = std::atan2(worldZ[1], worldZ[0]); 
    const int key = (ang + M_PI) / (2*M_PI) * kBinMaxCount;
    binV_[key] = bin;
  } else if (std::abs(localG[2]) > 0.9) {
    //  world vertical is approx. vertical
    const vec3 worldX = att.matrix() * vec3(1,0,0);
    const scalar_t ang = std::atan2(worldX[1], worldX[0]);
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
    vec3 bias(0,0,0), scl(1,1,1);

    std::vector<vec3> meas;
    for (const std::pair<int, SampleBin>& s : binH_) {
      meas.push_back(s.second.field);
    }
    for (const std::pair<int, SampleBin>&s : binV_) {
      meas.push_back(s.second.field);
    }
    const size_t N = meas.size();

    //  fit to sphere
    Matrix<scalar_t,Eigen::Dynamic,Eigen::Dynamic> A(N, 4);
    Matrix<scalar_t,Eigen::Dynamic,1> b(N, 1);
    for (size_t i=0; i < meas.size(); i++) {
      A(i,0) = 2*meas[i][0];
      A(i,1) = 2*meas[i][1];
      A(i,2) = 2*meas[i][2];
      A(i,3) = 1;
      b(i,0) = meas[i][0]*meas[i][0] + meas[i][1]*meas[i][1] + 
          meas[i][2]*meas[i][2];
    }
    //  solve system for center
    const Matrix<scalar_t,4,1> x = A.colPivHouseholderQr().solve(b);
    bias = x.block<3,1>(0,0);

    scalar_t mean_rad = 0.0, mean_rad_sqr = 0.0;
    //  calculate estimate of mean radius
    for (const vec3 &v : meas) {
      const scalar_t r = (v - bias).norm();
      mean_rad += r;
      mean_rad_sqr += r * r;
    }
    mean_rad /= N;
    mean_rad_sqr /= N;

    //  refine with GN-NLS
    Matrix<scalar_t, Eigen::Dynamic, 6> J(N, 6);
    Matrix<scalar_t, Eigen::Dynamic, 1> r(N, 1);
    Matrix<scalar_t, Eigen::Dynamic, Eigen::Dynamic> W;
    W.setZero();

    for (int iter = 0; iter < 20; iter++) {
      for (size_t i = 0; i < meas.size(); i++) {
        const scalar_t x = (meas[i][0] - bias[0]) / scl[0];
        const scalar_t y = (meas[i][1] - bias[1]) / scl[1];
        const scalar_t z = (meas[i][2] - bias[2]) / scl[2];
        const scalar_t r2 = x*x + y*y + z*z;
        
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
      scalar_t sigmaSquared=0;
      for (int j=0; j < r.rows(); j++) {
        sigmaSquared += r(j,0)*r(j,0);
      }
      sigmaSquared /= r.rows();
      //  calculate cauchy weights
      for (int j=0; j < r.rows(); j++) {
        const scalar_t errSqr = r(j,0)*r(j,0);
        W(j,j) = 1.0 / (1 + errSqr/sigmaSquared);
      }
      
      Matrix<scalar_t, 6, 6> H = J.transpose() * W * J;
      for (int i = 0; i < 3; i++) {  //  prior on scale
        H(i,i) *= 1.01;
      }
      Eigen::LDLT<Matrix<scalar_t,6,6>> LDLT(H);
      const Matrix<scalar_t, 6, 1> update = LDLT.solve(J.transpose() * W * r);

      scl += update.block<3, 1>(0, 0);
      bias += update.block<3, 1>(3, 0);
    }

    bias_ = bias;
    scale_ = scl;
  } else {
    bias_.setZero();
    scale_.setOnes();
  }
  
  calibrated_ = true;
}

} //  namespace kr
