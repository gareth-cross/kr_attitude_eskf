/*
 * AttitudeESKF.cpp
 *
 *  Copyright (c) 2013 Gareth Cross. All rights reserved.
 *
 *  This file is part of kr_attitude_eskf.
 *
 *	Created on: 12/24/2013
 *		  Author: gareth
 */

#define NDEBUG

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

#include "AttitudeESKF.hpp"

using namespace std;
using namespace Eigen;
using namespace kr;

//	skew symmetric matrix
template <typename T>
static inline Matrix<T, 3, 3> cross_skew(const Matrix<T, 3, 1> &w) {
  Matrix<T, 3, 3> W;

  W(0, 0) = 0;
  W(0, 1) = -w(2);
  W(0, 2) = w(1);

  W(1, 0) = w(2);
  W(1, 1) = 0;
  W(1, 2) = -w(0);

  W(2, 0) = -w(1);
  W(2, 1) = w(0);
  W(2, 2) = 0;

  return W;
}

//	hardcoded 3x3 invert (unchecked)
template <typename T>
static inline Matrix<T, 3, 3> invert(const Matrix<T, 3, 3> &A, T det) {
  Matrix<T, 3, 3> C;
  det = 1 / det;

  C(0, 0) = (-A(2, 1) * A(1, 2) + A(1, 1) * A(2, 2)) * det;
  C(0, 1) = (-A(0, 1) * A(2, 2) + A(0, 2) * A(2, 1)) * det;
  C(0, 2) = (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) * det;

  C(1, 0) = (A(2, 0) * A(1, 2) - A(1, 0) * A(2, 2)) * det;
  C(1, 1) = (-A(2, 0) * A(0, 2) + A(0, 0) * A(2, 2)) * det;
  C(1, 2) = (A(1, 0) * A(0, 2) - A(0, 0) * A(1, 2)) * det;

  C(2, 0) = (-A(2, 0) * A(1, 1) + A(1, 0) * A(2, 1)) * det;
  C(2, 1) = (A(2, 0) * A(0, 1) - A(0, 0) * A(2, 1)) * det;
  C(2, 2) = (-A(1, 0) * A(0, 1) + A(0, 0) * A(1, 1)) * det;

  return C;
}

//	hardcoded determinant
template <typename T> static inline T determinant(const Matrix<T, 3, 3> &A) {
  return A(0, 0) * (A(1, 1) * A(2, 2) - A(1, 2) * A(2, 1)) -
         A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) +
         A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));
}

AttitudeESKF::AttitudeESKF()
    : q_(), steadyCount_(0), biasThresh_(0), isStable_(true) {
  P_.setIdentity();
  b_.setZero();
  w_.setZero();
  dx_.setZero();

  magRef_.setZero();
  predMag_.setZero();

  estBias_ = false;
  useMag_ = false;
}

void AttitudeESKF::predict(const AttitudeESKF::vec3 &wb, AttitudeESKF::scalar_t dt) {
  static const Matrix<scalar_t, 3, 3> I3 =
      Matrix<scalar_t, 3, 3>::Identity(); //  identity R3

  scalar_t wb2 = wb[0]*wb[0] + wb[1]*wb[1] + wb[2]*wb[2];
  if (wb2 < biasThresh_*biasThresh_) {
    steadyCount_++; //  not rotating, update moving average

    if (estBias_ && steadyCount_ > 20) {
      b_ = (b_ * (steadyCount_ - 1) + wb) / steadyCount_;
    }
  } else {
    steadyCount_ = 0;
  }

  w_ = (wb - b_); //	true gyro reading

  //	error-state jacobian
  Matrix<scalar_t, 3, 3> F = I3 - cross_skew<scalar_t>(w_ * dt);
  
  //  integrate state and covariance
  q_.integrateEuler(quat<scalar_t>(0, w_[0], w_[1], w_[2]), dt);

  P_ = F * P_ * F.transpose();

  for (int i = 0; i < 3; i++) {
    P_(i, i) += var_.gyro[i];
  }
}

void AttitudeESKF::update(const AttitudeESKF::vec3 &ab, const AttitudeESKF::vec3 &mb) {
  Matrix<scalar_t, 3, 3> A;  //  for updating covariance

  //  rotation matrix: world -> body
  const Matrix<scalar_t, 3, 3> bRw = q_.conjugate().toMatrix();

  vec3 gravity;
  gravity[0] = 0.0;
  gravity[1] = 0.0;
  gravity[2] = 1.0;

  //  predicted gravity vector
  const vec3 aPred = bRw * gravity;

  if (!useMag_) {
    //  calculate jacobian
    Matrix<scalar_t, 3, 3> H = cross_skew(aPred);
    Matrix<scalar_t, 3, 1> r = ab - aPred;

    //  solve for the kalman gain
    Matrix<scalar_t, 3, 3> S = H * P_ * H.transpose();
    Matrix<scalar_t, 3, 3> Sinv;
    
    for (int i = 0; i < 3; i++) {
      S(i, i) += var_.accel[i];
    }

    const scalar_t det = determinant(S);
    if (std::abs(det) < static_cast<scalar_t>(1e-5)) {
      isStable_ = false;
      return;
    } else {
      isStable_ = true;
    }
    Sinv = invert(S, det);

    const Matrix<scalar_t, 3, 3> K = P_ * H.transpose() * Sinv;
    
    A = K * H;
    dx_ = K * r;
  }
  else {
#ifdef ATTITUDE_ESKF_BUILD_MAG  //  stop compilation of FullPivLU
    //  m-field prediction
    vec3 field = bRw * magRef_;
    predMag_ = field;

    Matrix<scalar_t, 6, 1> r;
    r.block<3, 1>(0, 0) = ab - aPred;
    r.block<3, 1>(3, 0) = mb - field;

    Matrix<scalar_t, 6, 3> H;
    H.setZero();

    //  jacobians for gravity and magnetic field
    H.block<3, 3>(0, 0) = cross_skew(aPred);
    H.block<3, 3>(3, 0) = cross_skew(field);

    //  covariance for both sensors
    Matrix<scalar_t, 6, 6> covR;
    covR.setZero();
    for (int i = 0; i < 3; i++) {
      covR(i, i) = var_.accel[i];
      covR(i + 3, i + 3) = var_.mag[i];
    }

    const Matrix<scalar_t, 6, 6> S = H * P_ * H.transpose() + covR;
    Matrix<scalar_t, 6, 6> Sinv;

    Eigen::FullPivLU<Matrix<scalar_t,6,6>> LU(S);
    isStable_ = LU.isInvertible();

    if (!isStable_) {
      return;
    }
    Sinv = LU.inverse();

    //  generate update
    const Matrix<scalar_t, 3, 6> K = P_ * H.transpose() * Sinv;
    dx_ = K * r;
    A = K * H;
#else
    dx_.setZero();
    A.setZero();
#endif
  }

  //  perform state update
  P_ = (Matrix<scalar_t, 3, 3>::Identity() - A) * P_;

  q_ = q_ * quat<scalar_t>(1.0, dx_[0], dx_[1], dx_[2]);
  q_ /= q_.norm();
}
