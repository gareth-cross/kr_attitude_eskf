/*
 * AttitudeESKF.cpp
 *
 *  Copyright (c) 2013 Gareth Cross. Apache 2 License.
 *
 *  This file is part of kr_attitude_eskf.
 *
 *	Created on: 12/24/2013
 *		  Author: gareth
 */

#ifndef NDEBUG
#define NDEBUG
#endif

#include "AttitudeESKF.hpp"
#include <Eigen/LU>
#include <Eigen/Cholesky>
#include <iostream>
#include <cmath>

using namespace Eigen;

namespace kr {

//	skew symmetric matrix
template <typename T>
static inline Matrix<T, 3, 3> crossSkew(const Matrix<T, 3, 1> &w) {
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

//  Eigen does not define these operators, which we use for integration
template <typename Scalar>
static inline Eigen::Quaternion<Scalar> operator + (const Eigen::Quaternion<Scalar>& a,
                                      const Eigen::Quaternion<Scalar>& b) {
  return Eigen::Quaternion<Scalar>(a.w()+b.w(),
                                   a.x()+b.x(),
                                   a.y()+b.y(),
                                   a.z()+b.z());
}

template <typename Scalar>
static inline Eigen::Quaternion<Scalar> operator * (const Eigen::Quaternion<Scalar>& q,
                                      Scalar s) {
  return Eigen::Quaternion<Scalar>(q.w() * s,
                                   q.x() * s,
                                   q.y() * s,
                                   q.z() * s);
}

/**
 *  @brief Integrate a rotation quaterion using Euler integration
 *  @param q Quaternion to integrate
 *  @param w Angular velocity (body frame), stored in 3 complex terms
 *  @param dt Time interval in seconds
 *  @param normalize If True, quaternion is normalized after integration
 */
template <typename Scalar>
static inline void integrateEuler(Eigen::Quaternion<Scalar> &q, Eigen::Quaternion<Scalar> &w, Scalar dt,
                    bool normalize = true) {
  q = q + (q * w * static_cast<Scalar>(0.5)) * dt;

  if (normalize) {
    q.normalize();
  }
}

/**
 *  @brief Integrate a rotation quaternion using 4th order Runge Kutta
 *  @param q Quaternion to integrate
 *  @param w Angular velocity (body frame), stored in 3 complex terms
 *  @param dt Time interval in seconds
 *  @param normalize If true, quaternion is normalized after integration
 */
template <typename Scalar>
static inline void integrateRungeKutta4(Eigen::Quaternion<Scalar> &q, const Eigen::Quaternion<Scalar> &w, Scalar dt,
                          bool normalize = true) {
  const static Scalar half = static_cast<Scalar>(0.5);
  const static Scalar two = static_cast<Scalar>(2);

  Eigen::Quaternion<Scalar> qw = q * w * half;
  Eigen::Quaternion<Scalar> k2 = (q + qw * dt * half) * w * half;
  Eigen::Quaternion<Scalar> k3 = (q + k2 * dt * half) * w * half;
  Eigen::Quaternion<Scalar> k4 = (q + k3 * dt) * w * half;

  q = q + (qw + k2 * two + k3 * two + k4) * (dt / 6);

  if (normalize) {
    q.normalize();
  }
}

template <typename Scalar>
static inline Eigen::Matrix<Scalar,3,3> 
rodrigues(const Eigen::Matrix<Scalar,3,1>& w) {
  const auto norm = w.norm();
  if (norm < std::numeric_limits<Scalar>::epsilon()*10) {
    return Eigen::Matrix<Scalar,3,3>::Identity() + crossSkew(w);
  }
  return Eigen::AngleAxis<Scalar>(norm, w / norm).matrix();
}

AttitudeESKF::AttitudeESKF()
    : q_(1,0,0,0), steadyCount_(0), biasThresh_(0), isStable_(true) {
  P_.setZero();
  b_.setZero();
  w_.setZero();
  dx_.setZero();

  magRef_.setZero();
  predMag_.setZero();

  estBias_ = false;
  ignoreZ_ = false;
  useMag_ = false;
}

void AttitudeESKF::predict(const AttitudeESKF::vec3 &wb,
                           AttitudeESKF::scalar_t dt, 
                           const AttitudeESKF::mat3 &cov,
                           bool useRK4) {
  static const Matrix<scalar_t, 3, 3> I3 =
      Matrix<scalar_t, 3, 3>::Identity(); //  identity R3

  scalar_t wb2 = wb[0] * wb[0] + wb[1] * wb[1] + wb[2] * wb[2];
  if (wb2 < biasThresh_ * biasThresh_) {
    steadyCount_++; //  not rotating, update moving average

    if (estBias_ && steadyCount_ > 20) {
      b_ = (b_ * (steadyCount_ - 1) + wb) / steadyCount_;
    }
  } else {
    steadyCount_ = 0;
  }

  w_ = (wb - b_); //	true gyro reading

  //	error-state jacobian
  const Matrix<scalar_t, 3, 3> F = I3 - crossSkew<scalar_t>(w_ * dt);

  //  integrate state and covariance
  Eigen::Quaternion<scalar_t> wQuat(0, w_[0], w_[1], w_[2]);
  if (!useRK4) {
    integrateEuler(q_, wQuat, dt, true);
  } else {
    integrateRungeKutta4(q_, wQuat, dt, true);
  }

  //  noise jacobian
  const Matrix <scalar_t,3,3> G = -I3 * dt;
  P_ = F*P_*F.transpose() + G*cov*G.transpose();
}

void AttitudeESKF::update(const AttitudeESKF::vec3 &ab, 
                          const mat3 &aCov, 
                          const AttitudeESKF::vec3 &mb, 
                          const mat3 &mCov) {
  Matrix<scalar_t, 3, 3> A;

  //  rotation matrix: world -> body
  const Matrix<scalar_t, 3, 3> bRw = q_.conjugate().matrix();

  vec3 gravity;
  gravity[0] = 0.0;
  gravity[1] = 0.0;
  gravity[2] = kOneG;

  //  predicted gravity vector
  const vec3 aPred = bRw * gravity;

  if (!useMag_) {
    //  calculate jacobian
    Matrix<scalar_t, 3, 3> H = crossSkew(aPred);
    Matrix<scalar_t, 3, 1> r = ab - aPred;

    //  solve for the kalman gain
    const Matrix<scalar_t, 3, 3> S = H * P_ * H.transpose() + aCov;
    Matrix<scalar_t, 3, 3> Sinv;

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
    H.block<3, 3>(0, 0) = crossSkew(aPred);
    H.block<3, 3>(3, 0) = crossSkew(field);

    //  covariance for both sensors
    Matrix<scalar_t, 6, 6> covR;
    covR.block<3,3>(0,0) = aCov;
    covR.block<3,3>(3,3) = mCov;

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
  
  if (ignoreZ_) {
    //  cancel body-frame z update
    dx_[2] = 0;
  }

  //  perform state update
  P_ = (Matrix<scalar_t, 3, 3>::Identity() - A) * P_;

  q_ = q_ * quat(1, dx_[0]/2, dx_[1]/2, dx_[2]/2);
  q_.normalize();
}
  
void AttitudeESKF::externalYawUpdate(scalar_t yaw, scalar_t alpha) {
  //  check if we are near the hover state
  const Matrix<scalar_t,3,3> wRb = q_.matrix();
  Matrix<scalar_t,3,1> g;
  g[0] = 0;
  g[1] = 0;
  g[2] = 1;
  
  g = wRb.transpose() * g;
  if (g[2] > 0.85) {
    //  break into roll pitch yaw
    Matrix<scalar_t,3,1> rpy = getRPY(wRb);
    //  interpolate between prediction and estimate
    rpy[2] = rpy[2]*(1-alpha) + yaw*alpha;
    q_ = Eigen::AngleAxis<scalar_t>(rpy[2],vec3(0,0,1)) *
    Eigen::AngleAxis<scalar_t>(rpy[1],vec3(0,1,0)) *
    Eigen::AngleAxis<scalar_t>(rpy[0],vec3(1,0,0));
  }
}

bool AttitudeESKF::initialize(const vec3 &ab,
                              const vec3 &aCov,
                              const vec3 &mb,
                              const vec3 &mCov) {
  if (!useMag_) {
    //  determine attitude angles
    scalar_t ay = ab[1];
    if (ay > kOneG) { ay = kOneG; }
    else if (ay < -kOneG) { ay = -kOneG; }
    const scalar_t& ax = ab[0];
    const scalar_t& az = ab[2]; 
    
    const scalar_t phi = std::asin(-ay / kOneG);  //  roll
    const scalar_t theta = std::atan2(ax, az);    //  pitch
  
    q_ = Eigen::AngleAxis<scalar_t>(theta, vec3(0,1,0)) * 
         Eigen::AngleAxis<scalar_t>(phi, vec3(1,0,0));
  }
  else {
    ///  @todo: This is kind of ugly, find some simpler mechanism to do this.
    
#ifdef ATTITUDE_ESKF_BUILD_MAG
    const static scalar_t eps(1e-6);
    for (int i=0; i < 3; i++) {
      if (aCov[i] < eps || mCov[i] < eps) {
        return false;
      }
    }
    //  jacobian
    Eigen::Matrix <scalar_t,6,3> J;
    J.block<3,3>(0,0) = crossSkew(ab);
    J.block<3,3>(3,0) = crossSkew(mb);
    
    //  weight matrix
    Eigen::Matrix <scalar_t,6,6> S;
    S.setZero();
    for (int i=0; i < 3; i++) {
      S(i,i) = 1 / aCov[i];
      S(i+3,i+3) = 1 / mCov[i];
    }
    
    //  hessian
    const mat3 H = J.transpose() * S * J;
    const Eigen::LDLT<mat3> ldlt(H);
    
    //  optimize
    vec3 w(0,0,0);
    Matrix<scalar_t,6,1> r;
    for (unsigned int iter=0; iter < 5; iter++) {
      const mat3 W = rodrigues(w);
      //  residuals
      r.block<3,1>(0,0) = (W * vec3(0,0,kOneG)) - ab;
      r.block<3,1>(3,0) = (W * magRef_) - mb;
      //  step
      w.noalias() += ldlt.solve(J.transpose() * S * r);
    }
    q_ = quat(rodrigues(w).transpose());
#endif
  }
  //  start w/ a large uncertainty
  P_.setIdentity();
  P_ *= M_PI*M_PI;
  
  return true;
}
  
AttitudeESKF::vec3 AttitudeESKF::getRPY(const mat3& R) {
  vec3 rpy;
  scalar_t sth = -R(2, 0);
  if (sth > 1) {
    sth = 1;
  } else if (sth < -1) {
    sth = -1;
  }
  
  const scalar_t theta = std::asin(sth);
  const scalar_t cth = std::sqrt(1 - sth*sth);
  
  scalar_t phi, psi;
  if (cth < static_cast<scalar_t>(1.0e-6)) {
    phi = std::atan2(R(0, 1), R(1, 1));
    psi = 0;
  } else {
    phi = std::atan2(R(2, 1), R(2, 2));
    psi = std::atan2(R(1, 0), R(0, 0));
  }
  
  rpy[0] = phi;    //  x, [-pi,pi]
  rpy[1] = theta;  //  y, [-pi/2,pi/2]
  rpy[2] = psi;    //  z, [-pi,pi]
  return rpy;
}

} //  namespace kr

