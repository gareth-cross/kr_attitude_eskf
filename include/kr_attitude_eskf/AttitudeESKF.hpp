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

#include <Eigen/Core>
#include <Eigen/Dense>

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

  typedef Eigen::Matrix<scalar_t, 3, 1> vec3; /// Vector in R3
  typedef Eigen::Matrix<scalar_t, 3, 3> mat3; /// Matrix in R3
  typedef Eigen::Quaternion<scalar_t> quat;   /// Member of S4

  static constexpr scalar_t kOneG = 9.80665;  /// Earth gravity
  
  /**
   *  @brief Ctor, initializes state to all zeros.
   */
  AttitudeESKF();

  /**
   * @brief predict Perform the prediction step.
   * @param wg Uncorrected gyroscope readings in body frame.
   * @param dt Time step in seconds.
   * @param cov Covariance on gyro measurement, units of (rad/s)^2.
   * @param useRK4 If true, use RK4 integration - otherwise euler is used.
   */
  void predict(const vec3 &wg, scalar_t dt, 
               const mat3& cov, bool useRK4=false);

  /**
   * @brief update Perform the update step.
   * @param ab Accelerometer reading in body frame (units of m/s^2).
   * @param aCov Covariance on accelerometer measurement, units of (m/s^2)^2.
   * @param mb Measured magnetic field in body frame (units of gauss)
   * @param mCov Covariance on magnetometer measurement, units of gauss^2.
   * 
   * @note Magnetometer elements only used if usesMagnetometer is set to true.
   */
  void update(const vec3 &ab,
              const mat3 &aCov,
              const vec3 &mb = vec3::Zero(),
              const mat3 &mCov = mat3::Zero());
  
  /**
   * @brief initialize Initialize the pose.
   * 
   * @param ab Measured body acceleration (units of m/s^2).
   * @param aCov Diagonal uncertainty on acceleration, units of (m/s^2)^2.
   * @param mb Measured magnetic field in body frame (units of gauss).
   * @param mCov Diagonal uncertainty on magnetic field, units of gauss^2.
   * @param maxIterations Max iterations during initialization.
   * 
   * @note Uses non-linear least squares to formulate initial rotation vector.
   * The inverse covariance matrices are used to weight the input vectors.
   * 
   * @note If magnetometer is disabled, a roll and pitch angles are determined
   * from the gravity vector. The yaw angle is zeroed.
   * 
   * @return True on success, false on failure. Failure will occur if aCov or
   * mCov are not invertible (and the magnetometer is enabled).
   */
  bool initialize(const vec3 &ab,
                  const vec3 &aCov,
                  const vec3 &mb = vec3::Zero(),
                  const vec3 &mCov = vec3::Zero(),
                  unsigned int maxIterations = 5);
  
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
   * @brief setMagneticReference Set the magnetic reference vector.
   * @param magRef Vector determining the North magnetic field.
   */
  void setMagneticReference(const vec3 &magRef) { magRef_ = magRef; }

  /**
   * @brief getQuat Get the state as a quaternion.
   * @return Instance of kr::quat.
   */
  const quat &getQuat() const { return q_; }

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
  const mat3 &getCovariance() const { return P_; }
  mat3 &getCovariance() { return P_; }
  
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
    
  quat q_; /// Orientation
  mat3 P_; /// System covariance

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
};

} // namespace kr

#endif /* defined(KR_ATTITUDE_ESKF_H_) */
