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

#include <kr_attitude_eskf/AttitudeESKF.hpp>
#include <kr_attitude_eskf/AttitudeMagCalib.hpp>

namespace kr_attitude_eskf {

class Node {
public:
  Node(const ros::NodeHandle& nh);
  
private:
  ros::NodeHandle nh_;
};

} //  namespace kr_attitude_eskf

#endif // KR_ATTITUDE_NODE_HPP_
