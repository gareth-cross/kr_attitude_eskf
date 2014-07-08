/*
 * AttitudeESKF.hpp
 *
 *  Copyright (c) 2013 Gareth Cross. All rights reserved.
 *
 *  This file is part of kr_attitude_eskf.
 *
 *	Created on: 12/24/2013
 *		  Author: gareth
 */

#ifndef KR_ATTITUDE_ESKF_H_
#define KR_ATTITUDE_ESKF_H_

#include <kr_math/quaternion.hpp>

namespace kr {

/**
 *  @class AttitudeESKF
 *  @brief Implementation of an error-state EKF for attitude determination using
 * quaternions.
 *  @note Gravity and magnetic field are the supported reference vectors.
 *  @see 'Attitude Error Representations for Kalman Filtering', F. Landis
 * Markley, JPL
 */
class AttitudeESKF {
public:
  typedef double scalar_t; /**< Type used for all calculations, change as
                              performance requires */

  typedef Eigen::Matrix<scalar_t, 3, 1> vec3; /**< Vector in R3 */

  /**
   * @brief Describes all the sensor noise properties.
   */
  struct VarSettings {
    scalar_t accel[3]; /// XYZ variance on acceleration, units of Gs
    scalar_t gyro[3];  /// XYZ variance on gyroscope, units of rad/s
    scalar_t mag[3];   /// XYZ variance on magnetometer, units of Gauss

    VarSettings() {
      for (int i = 0; i < 3; i++) {
        accel[i] = gyro[i] = mag[i] = 0;
      }
    }
  };

  /**
   *  @brief Ctor, initializes state to all zeros.
   */
  AttitudeESKF();

  /**
   *  @brief predict Perform the prediction step.
   *  @param wg Uncorrected gyroscope readings in body frame.
   *  @param dt Time step in seconds.
   *
   *  @note Integrates the nominal state using RK4.
   */
  void predict(const vec3 &wg, scalar_t dt);

  /**
   *  @brief update Perform the update step.
   *  @param ab Accelerometer reading in body frame (units of Gs).
   */
  void update(const vec3 &ab, const vec3 &mb = vec3::Zero());

  /**
   * @brief setEstimatesBias Enable/Disable bias estimation.
   * @param estBias If true, bias is estimated online.
   * @note Biases are estimated with a moving average filter when the platform is not in motion.
   */
  void setEstimatesBias(bool estBias) { estBias_ = estBias; }

  /**
   * @brief setGyroBiasThreshold Set the gyro bias threshold.
   * @param thresh Threshold, units of rad/s.
   * @note Below this threshold of rotation, gyro biases are estimated.
   */
  void setGyroBiasThreshold(scalar_t thresh) { biasThresh_ = thresh; }

  /**
   * @brief setUsesMagnetometer Enable/disable magnetometer update.
   * @param useMag If true, mag update is enabled.
   */
  void setUsesMagnetometer(bool useMag) { useMag_ = useMag; }

  /**
   * @brief setVariances Set the process/measurement variances.
   * @param var Instance of VarSettings.
   */
  void setVariances(const VarSettings &var) { var_ = var; }

  /**
   * @brief setMagneticReference Set the magnetic reference vector.
   * @param magRef Vector determining the North magnetic field.
   */
  void setMagneticReference(const vec3 &magRef) { magRef_ = magRef; }

  /**
   * @brief getQuat Get the state as a quaternion.
   * @return Instance of kr::quat.
   */
  const kr::quat<scalar_t> &getQuat() const { return q_; }

  /**
   * @brief getAngularVelocity Get angular velocity (corrected for bias).
   * @return Angular velocity in rad/s.
   */
  const vec3 &getAngularVelocity() const { return w_; }

  /**
   * @brief getGyroBias Get the current gyro bias estimate.
   * @return Gyro bias in units of rad/s.
   */
  const vec3 &getGyroBias() const { return b_; }

  /**
   * @brief getCovariance Get the system covariance on the error state.
   * @return 3x3 covariance matrix.
   */
  const Eigen::Matrix<scalar_t, 3, 3> &getCovariance() const { return P_; }

  /**
   * @brief getPredictedField Get the predicted magnetic field.
   * @return The predicted magnetic field for the current state, units of gauss.
   */
  const vec3 &getPredictedField() const { return predMag_; }

  /**
   *  @brief getCorrection Get the last correction (error state) generated.
   *  @return The previous error state.
   */
  const vec3 &getCorrection() const { return dx_; }

  /**
   * @brief isStable Determine if the filter is stable.
   * @return True if the Kalman gain was non-singular at the last update.
   */
  bool isStable() const { return isStable_; }

private:
  kr::quat<scalar_t> q_;            /// Orientation
  Eigen::Matrix<scalar_t, 3, 3> P_; /// System covariance

  vec3 w_;
  vec3 b_;
  unsigned long steadyCount_;
  scalar_t biasThresh_;

  vec3 magRef_;
  vec3 predMag_;

  vec3 dx_;

  bool isStable_;
  bool estBias_;
  bool useMag_;

  VarSettings var_;
};

} // namespace kr

#endif /* defined(KR_ATTITUDE_ESKF_H_) */
