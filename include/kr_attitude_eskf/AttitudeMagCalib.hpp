/*
 * AttitudeMagCalib.hpp
 *
 *  Copyright (c) 2013 Gareth Cross. All rights reserved.
 *
 *  This file is part of kr_attitude_eskf.
 *
 *	Created on: 23/06/2014
 *		  Author: gareth
 */

#ifndef KR_ATTITUDE_MSG_CALIB_H_
#define KR_ATTITUDE_MSG_CALIB_H_

#include <map>
#include <vector>
#include "kr_attitude_eskf/AttitudeESKF.hpp"

namespace kr {

/**
 * @class AttitudeMagCalib
 * @brief Tool for calibrating the system magnetometer.
 */
class AttitudeMagCalib {
public:
  /// Types
  typedef kr::AttitudeESKF::scalar_t scalar_t;
  typedef kr::AttitudeESKF::vec3 vec3;
  typedef kr::AttitudeESKF::mat3 mat3;
  typedef kr::AttitudeESKF::quat quat;
  
  /**
   * @brief Types of calibration that can be performed.
   */
  enum CalibrationType : int {
    FullCalibration = 0,      /// Calculate bias and scale
  };

  /**
   * @brief singular_hessian Thrown when non-linear least squares fails because the hessian cannot be inverted.
   */
  class singular_hessian : public std::exception {};

  /**
   * @brief insufficient_data Thrown if calibration is attempted with an insufficient number of points.
   */
  class insufficient_data : public std::exception {};

  /**
   * @brief Ctor, initializes to reset state.
   */
  AttitudeMagCalib();

  /**
   * @brief reset Erase collected points and reset the calibration state.
   */
  void reset();

  /**
   * @brief appendSample Append an attitude estimate and field measurement.
   * @param att Current attitude at time of measurement.
   * @param field Current magnetic field.
   */
  void appendSample(const quat &att, const vec3 &field);

  /**
   * @brief isReady Is the system ready to calibrate?
   * @return True if enough points have been collected.
   */
  bool isReady() const;

  /**
   * @brief isCalibrated Are the calibration values set?
   * @return True if calibrate() ran succesfully and the system has not been reset.
   */
  bool isCalibrated() const;

  /**
   * @brief Run calibration on data recorded.
   * @throws insufficient_data, singular_hessian
   * 
   * @note This method will optimize the scale and bias parameters using
   * a two step process:
   *  1) Least-squares fit of a sphere as an initial guess.
   *  2) LM fit of an axis-aligned spheroid as a refined estimate.
   */
  void calibrate(AttitudeMagCalib::CalibrationType type = FullCalibration);

  /**
   * @brief getBias Get the estimate of magnetometer bias.
   * @return Vector in R3, units of gauss.
   */
  const vec3& getBias() const { return bias_; }

  /**
   * @brief getScale Get the estimate of magnetometer scale.
   * @return Vector in R3, unitless.
   */
  const vec3& getScale() const { return scale_; }

private:
  constexpr static int kBinMaxCount = 40;   /// Number of samples to take on each axis

  struct SampleBin {
    vec3 field; ///  Field measured in this sample
    quat q;     ///  Unreferenced quaternion for this sample
  };

  std::map<int, SampleBin> binV_;
  std::map<int, SampleBin> binH_;
  bool calibrated_;

  vec3 bias_;
  vec3 scale_;
};

} // namespace kr

#endif // KR_ATTITUDE_MSG_CALIB_H_
