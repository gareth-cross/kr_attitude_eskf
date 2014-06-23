/*
 * AttitudeESKF.cpp
 *
 *  Copyright (c) 2013 Gareth Cross. All rights reserved.
 *
 *  This file is part of AttitudeESKF.
 *
 *	Created on: 12/24/2013
 *		  Author: gareth
 */

#include "AttitudeESKF.hpp"
#include "SO3.hpp"

using namespace std;
using namespace Eigen;
using namespace kr;

//	skew symmetric matrix
template <typename T>
static inline Matrix<T,3,3> cross_skew(const Matrix<T,3,1>& w)
{
    Matrix<T,3,3> W = Matrix<T,3,3>::Zero();

    W(0,1) = -w(2);
    W(0,2) = w(1);

    W(1,0) = w(2);
    W(1,2) = -w(0);

    W(2,0) = -w(1);
    W(2,1) = w(0);

    return W;
}

//	hardcoded 3x3 invert (unchecked)
template <typename T>
static inline Matrix<T,3,3> invert(const Matrix<T,3,3>& A, T det)
{
    Matrix<T,3,3> C;
    det = 1.0 / det;

    C(0,0) = (-A(2,1)*A(1,2) + A(1,1)*A(2,2)) * det;
    C(0,1) = (-A(0,1)*A(2,2) + A(0,2)*A(2,1)) * det;
    C(0,2) = ( A(0,1)*A(1,2) - A(0,2)*A(1,1)) * det;

    C(1,0) = ( A(2,0)*A(1,2) - A(1,0)*A(2,2)) * det;
    C(1,1) = (-A(2,0)*A(0,2) + A(0,0)*A(2,2)) * det;
    C(1,2) = ( A(1,0)*A(0,2) - A(0,0)*A(1,2)) * det;

    C(2,0) = (-A(2,0)*A(1,1) + A(1,0)*A(2,1)) * det;
    C(2,1) = ( A(2,0)*A(0,1) - A(0,0)*A(2,1)) * det;
    C(2,2) = (-A(1,0)*A(0,1) + A(0,0)*A(1,1)) * det;

    return C;
}

//	hardcoded determinant
template <typename T>
static inline T determinant(const Matrix<T,3,3>& A)
{
    return  A(0,0) * ( A(1,1)*A(2,2) - A(1,2)*A(2,1) ) -
            A(0,1) * ( A(1,0)*A(2,2) - A(1,2)*A(2,0) ) +
            A(0,2) * ( A(1,0)*A(2,1) - A(1,1)*A(2,0) );
}

AttitudeESKF::AttitudeESKF() :
  q_(), isStable_(true), lastTime_(0.0), steadyCount_(0), biasThresh_(0.0)
{
  P_.setZero();
  b_.setZero();
  w_.setZero();

  magRef_.setZero();
  predMag_.setZero();

  estBias_ = false;
  useMag_ = false;
}

void AttitudeESKF::predict(const vec3 &wb, double time)
{
  static const Matrix<scalar_t,3,3> I3 = Matrix<scalar_t,3,3>::Identity(); //  identity R3

  scalar_t dt = 0.01;
  if (lastTime_ != 0.0 && time > lastTime_) {
    dt = static_cast<scalar_t>(time - lastTime_);
  }
  lastTime_ = time;

  if (wb.norm() < biasThresh_) {
    steadyCount_++; //  not rotating, update moving average

    if (estBias_ && steadyCount_ > 20)
    {
      b_ = (b_ * (steadyCount_-1) + wb) / steadyCount_;
    }
  } else {
    steadyCount_ = 0;
  }

  w_ = (wb - b_);   //	true gyro reading

  //	error-state jacobian
  Matrix<scalar_t,3,3> F;
  F.setZero();
  F.block<3,3>(0,0) = I3 - cross_skew<scalar_t>(w_ * dt);

  //  integrate nominal state
  q_.integrateRungeKutta4(quat<scalar_t>(0.0, w_[0], w_[1], w_[2]), dt);

  //  integrate covariance
  P_ = F * P_ * F.transpose();

  for (int i=0; i < 3; i++) {
    P_(i,i) += var_.gyro[i];
  }
}

void AttitudeESKF::update(const vec3 &ab, const vec3 &mb)
{
  Matrix<scalar_t,3,1> dx;  //  error state
  Matrix<scalar_t,3,3> A;   //  for updating covariance

  //  rotation matrix: body -> world
  const Matrix<scalar_t,3,3> R = q_.to_matrix().cast<scalar_t>();

  vec3 gravity;
  gravity[0] = 0.0;
  gravity[1] = 0.0;
  gravity[2] = 1.0;

  //  predicted gravity vector
  const vec3 aPred = R.transpose() * gravity;

  if (!useMag_)
  {
    Matrix <scalar_t,3,3> covR;
    covR.setZero();
    for (int i=0; i < 3; i++) {
      covR(i,i) = var_.accel[i];
    }

    //  calculate jacobian
    Matrix <scalar_t,3,3> H = cross_skew(aPred);
    Matrix<scalar_t,3,1> r = ab - aPred;

    //  solve for the kalman gain
    Matrix<scalar_t,3,3> S = H * P_ * H.transpose() + covR;
    Matrix<scalar_t,3,3> Sinv;

    const scalar_t det = determinant(S);
    if (std::abs(det) < 1e-5) {
      isStable_ = false;
      return;
    } else {
      isStable_ = true;
    }
    Sinv = invert(S,det);

    const Matrix<scalar_t,3,3> K = P_ * H.transpose() * Sinv;

    dx = K * r;
    A = K * H;
  }
  else
  {
    //  m-field prediction
    vec3 field = R.transpose() * magRef_;
    predMag_ = field;

    Matrix<scalar_t,6,1> r;
    r.block<3,1>(0,0) = ab - aPred;
    r.block<3,1>(3,0) = mb - field;

    Matrix <scalar_t,6,3> H;
    H.setZero();

    //  jacobians for gravity and magnetic field
    H.block<3,3>(0,0) = cross_skew(aPred);
    H.block<3,3>(3,0) = cross_skew(field);

    //  covariance for both sensors
    Matrix <scalar_t,6,6> covR;
    covR.setZero();
    for (int i=0; i < 3; i++) {
      covR(i,i) = var_.accel[i];
      covR(i+3,i+3) = var_.mag[i];
    }

    const Matrix<scalar_t,6,6> S = H * P_ * H.transpose() + covR;
    Matrix<scalar_t,6,6> Sinv;

    auto LU = S.fullPivLu();
    isStable_ = LU.isInvertible();

    if (!isStable_) {
        return;
    }
    Sinv = LU.inverse();

    //  generate update
    const Matrix<scalar_t,3,6> K = P_ * H.transpose() * Sinv;
    dx = K * r;
    A = K * H;
  }

  //  perform state update
  P_ = (Matrix<scalar_t,3,3>::Identity() - A) * P_;

  q_ = q_ * quat<scalar_t>(1.0, dx[0], dx[1], dx[2]);
  q_ /= q_.norm();
}


//  (world) = Rz * Ry * Rx (body)
Eigen::Matrix<AttitudeESKF::scalar_t,3,1> AttitudeESKF::getRPY() const
{
	const Matrix<scalar_t,3,3> R = q_.to_matrix();
	Matrix<scalar_t,3,1> rpy;

	scalar_t sth = -R(2,0);
	if (sth > 1.0) {
		sth = 1.0;
	} else if (sth < -1.0) {
		sth = -1.0;
	}

	const scalar_t theta = std::asin(sth);
	const scalar_t cth = std::sqrt(1.0 - sth*sth);

	scalar_t phi, psi;
	if (cth < 1e-6)
	{
		phi = std::atan2(R(0,1), R(1,1));
		psi = 0.0;
	}
	else
	{
		phi = std::atan2(R(2,1), R(2,2));
		psi = std::atan2(R(1,0), R(0,0));
	}

	rpy[0] = phi;   //  x
	rpy[1] = theta; //  y
	rpy[2] = psi;   //  z
	return rpy;
}
