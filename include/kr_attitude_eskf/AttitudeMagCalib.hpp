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
#include <kr_math/base_types.hpp>

namespace kr {

/**
 *  @class AttitudeMagCalib
 *  @brief Tool for calibrating the system magnetometer.
 */
class AttitudeMagCalib {
public:
  /**
   * @brief Types of calibration that can be performed.
   */
  enum CalibrationType : int {
    FullCalibration = 0,      /// Calculate bias, scale and reference
    ReferenceCalibration = 1, /// Calibrate only North reference vector
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
  void appendSample(const kr::quatd &att, const Eigen::Vector3d &field);

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
   * @brief calibrate
   * @throws
   */
  void calibrate(AttitudeMagCalib::CalibrationType type);

  /**
   * @brief getBias Get the estimate of magnetometer bias.
   * @return Vector in R3, units of gauss.
   */
  const Eigen::Vector3d& getBias() const { return bias_; }

  /**
   * @brief getScale Get the estimate of magnetometer scale.
   * @return Vector in R3, unitless.
   */
  const Eigen::Vector3d& getScale() const { return scale_; }

  /**
   * @brief getReference Get the magnetic north vector.
   * @return Vector in R3, units of gauss.
   *
   * @note The north reference vector corresponds to the vector [H 0 V], where H is the horizontal component of the B-field, and V the vertical component.
   */
  const Eigen::Vector3d& getReference() const { return ref_; }

private:
  constexpr static int kBinMaxCount = 30;   /// Number of samples to take on each axis

  struct SampleBin {
    Eigen::Vector3d field; ///  Field measured in this sample
    kr::quatd q;           ///  Unreferenced quaternion for this sample
  };

  std::map<int, SampleBin> binV_;
  std::map<int, SampleBin> binH_;
  bool calibrated_;

  Eigen::Vector3d bias_;
  Eigen::Vector3d scale_;
  Eigen::Vector3d ref_;
};

} // namespace kr

#endif // KR_ATTITUDE_MSG_CALIB_H_
