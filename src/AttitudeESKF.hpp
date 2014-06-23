/*
 * AttitudeESKF.hpp
 *
 *  Copyright (c) 2013 Gareth Cross. All rights reserved.
 *
 *  This file is part of AttitudeESKF.
 *
 *	Created on: 12/24/2013
 *		  Author: gareth
 */

#ifndef __AttitudeESKF__
#define __AttitudeESKF__

#include "quaternion.hpp"

namespace kr {

/**
 *  @class AttitudeESKF
 *  @brief Implementation of an error-state EKF for attitude determination using Quaternions.
 *  @note Gravity and magnetic field are the supported reference vectors.
 *  @see 'Attitude Error Representations for Kalman Filtering', F. Landis Markley, JPL
 */
class AttitudeESKF
{
public:

  typedef double scalar_t;  /**< Type used for all calculations, change as performance requires */

  typedef Eigen::Matrix<scalar_t,3,1> vec3; /**< Vector in R3 */

  /**
   * @brief The VarSettings struct
   */
  struct VarSettings {
    scalar_t accel[3];  /// XYZ variance on acceleration, units of Gs
    scalar_t gyro[3];   /// XYZ variance on gyroscope, units of rad/s
    scalar_t mag[3];    /// XYZ variance on magnetometer, units of Gauss

    VarSettings() {
      for (int i=0; i < 3; i++) {
        accel[i] = gyro[i] = mag[i] = 0.0;
      }
    }
  };

  /**
   *  @brief Ctor, initializes state to all zeros
   */
  AttitudeESKF();

  /**
   *  @brief Perform the prediction step
   *  @param wg Uncorrected gyroscope readings in body frame
   *  @param time Current time in seconds
   *
   *  @note Integrates the nominal state using RK4.
   */
  void predict(const vec3& wg, double time);

  /**
   *  @brief Perform the update step
   *  @param ab Accelerometer reading in body frame (units of Gs)
   */
  void update(const vec3& ab, const vec3& mb = vec3::Zero());

  /**
   *	@brief Get Roll-Pitch-Yaw as a 3-element vector
   */
  Eigen::Matrix<scalar_t,3,1> getRPY() const;

  /**
   * @brief setEstimatesBias
   * @param estBias
   */
  void setEstimatesBias(bool estBias) { estBias_ = estBias; }

  /**
   * @brief setGyroBiasThreshold
   * @param thresh
   */
  void setGyroBiasThreshold(scalar_t thresh) { biasThresh_ = thresh; }

  /**
   * @brief setUsesMagnetometer
   * @param useMag
   */
  void setUsesMagnetometer(bool useMag) { useMag_ = useMag; }

  /**
   * @brief setVariances
   * @param var
   */
  void setVariances(const VarSettings& var) { var_ = var; }

  /**
   * @brief setMagneticReference
   * @param magRef
   */
  void setMagneticReference(const vec3& magRef) { magRef_ = magRef; }

  /**
   * @brief getQuat
   * @return
   */
  const kr::quat<scalar_t>& getQuat() const { return q_; }

  /**
   * @brief getAngularVelocity
   * @return
   */
  const vec3& getAngularVelocity() const { return w_; }

  /**
   * @brief getGyroBias
   * @return
   */
  const vec3& getGyroBias() const { return b_; }

  /**
   * @brief getCovariance
   * @return
   */
  const Eigen::Matrix<scalar_t,3,3>& getCovariance() const { return P_; }

  /**
   * @brief
   * @return
   */
  const vec3& getPredictedField() const { return predMag_; }

  /**
   * @brief isStable
   * @return
   */
  bool isStable() const { return isStable_; }

private:
  kr::quat<scalar_t> q_;           /// Orientation
  Eigen::Matrix<scalar_t,3,3> P_;  /// System covariance
  double lastTime_;

  vec3 w_;
  vec3 b_;
  unsigned long steadyCount_;
  scalar_t biasThresh_;

  vec3 magRef_;  //  North
  vec3 predMag_;

  bool isStable_;
  bool estBias_;
  bool useMag_;

  VarSettings var_;
};

}  // namespace kr

#endif /* defined(__AttitudeESKF__) */
