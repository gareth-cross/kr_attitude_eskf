/*
 * Node.hpp
 *
 *  Copyright (c) 2014 Gareth Cross. All rights reserved.
 *
 *  This file is part of kr_attitude_eskf.
 *
 *	Created on: 31/08/2014
 *		  Author: gareth
 */

#ifndef KR_ATTITUDE_NODE_HPP_
#define KR_ATTITUDE_NODE_HPP_

#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/MagneticField.h>
#include <std_srvs/Empty.h>
#include <diagnostic_updater/diagnostic_updater.h>
#include <diagnostic_updater/publisher.h>
#include <message_filters/synchronizer.h>
#include <message_filters/subscriber.h>
#include <message_filters/sync_policies/exact_time.h>

#include <kr_math/base_types.hpp>
#include <kr_attitude_eskf/AttitudeESKF.hpp>
#include <kr_attitude_eskf/AttitudeMagCalib.hpp>

namespace kr_attitude_eskf {

class Node {
public:
  Node(const ros::NodeHandle& nh, const ros::NodeHandle& pnh);
  
private:
  static constexpr unsigned int kROSQueueSize = 100;
  
  ros::NodeHandle nh_;
  ros::Publisher pubImu_;
  ros::Publisher pubBias_;
  ros::Publisher pubPose_;
  ros::ServiceServer srvCalibrate_;
  
  //  subsribers
  message_filters::Subscriber<sensor_msgs::Imu> subImu_;
  message_filters::Subscriber<sensor_msgs::MagneticField> subField_;
  ros::Subscriber subImuUnsync_;
  
  //  time sync
  typedef message_filters::sync_policies::ExactTime<sensor_msgs::Imu,
    sensor_msgs::MagneticField> TimeSyncPolicy;
  message_filters::Synchronizer<TimeSyncPolicy> sync_;
  
  //  options
  bool enableMag_;
  double gyroBiasThresh_;
  
  //  implementation
  kr::vec3d magBias_;
  kr::vec3d magScale_;
  kr::vec3d magReference_;
  kr::AttitudeESKF eskf_;
  kr::AttitudeMagCalib calib_;
  ros::Time prevStamp_;
  enum {
    Uncalibrated = 0,
    Calibrating = 1,
    CalibrationComplete = 2,
  } calibState_;
  
  //  callbacks
  void inputCallback(const sensor_msgs::ImuConstPtr&,
                     const sensor_msgs::MagneticFieldConstPtr&);
  
  bool beginCalibration(std_srvs::Empty::Request&,
                        std_srvs::Empty::Response&);
};

} //  namespace kr_attitude_eskf

#endif // KR_ATTITUDE_NODE_HPP_
