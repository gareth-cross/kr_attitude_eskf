/*
 * kr_attitude_eskf.cpp
 *
 *  Copyright (c) 2014 Kumar Robotics. All rights reserved.
 *
 *  This file is part of kr_attitude_eskf.
 *
 *  Created on: 17/6/2014
 *		  Author: gareth
 */

#include <kr_attitude_eskf/Node.hpp>

int main(int argc, char **argv) {
  ros::init(argc, argv, "kr_attitude_eskf");
  ros::NodeHandle nh;
  ros::NodeHandle pnh("~");
  kr_attitude_eskf::Node node(nh,pnh);
  ros::spin();
  return 0;
}
