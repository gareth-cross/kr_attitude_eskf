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

#include <ros/ros.h>
#include <sensor_msgs/Imu.h>
#include <sensor_msgs/MagneticField.h>
#include <sensor_msgs/FluidPressure.h>

#include <geometry_msgs/Vector3Stamped.h>
#include <geometry_msgs/PoseStamped.h>

#include <eigen_conversions/eigen_msg.h>

#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>

#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>

#include <kr_attitude_eskf/Status.h>

#include "AttitudeESKF.hpp"
#include "AttitudeMagCalib.hpp"

#include <kr_math/SO3.hpp>

#include <algorithm>
#include <string>

using namespace std;
using namespace Eigen;

kr::AttitudeESKF eskf;
kr::AttitudeMagCalib calib;
kr::AttitudeMagCalib::CalibrationType preferredCalib =
    kr::AttitudeMagCalib::FullCalibration;

ros::Publisher pubImu;
ros::Publisher pubBias;
ros::Publisher pubField;
ros::Publisher pubStatus;
ros::Publisher pubPose;

ros::NodeHandlePtr nh; //  main node handle

const std::string sfx[] = { "x", "y", "z" }; //  suffixes to parameter names

boost::shared_ptr<tf::TransformBroadcaster> tfBroadcaster;
bool broadcast_frame = false;

std::string bodyFrameName;  //  name of frame_id
bool publish_pose = false;

enum {
  MagIdle = 0,
  MagInCalibration = 1,
  MagCalibrated = 2,
} topicMode = MagIdle;

bool subscribe_mag = false; //  use the magnetometer or not?

//  default bias and scale
Vector3d magBias = Vector3d::Zero();
Vector3d magScale = Vector3d::Ones();
Vector3d magReference = Vector3d::Zero();

//  last timestamp
ros::Time prevStamp(0,0);

void imu_callback(const sensor_msgs::ImuConstPtr &imu,
                  const sensor_msgs::MagneticFieldConstPtr &field) {

  Vector3d wm; //  measured angular rate
  Vector3d am; //  measured acceleration
  Vector3d mm; //  measured magnetic field

  tf::vectorMsgToEigen(imu->angular_velocity, wm);
  tf::vectorMsgToEigen(imu->linear_acceleration, am);

  if (subscribe_mag) {
    tf::vectorMsgToEigen(field->magnetic_field, mm);
  } else {
    mm.setZero(); //  safe default
  }

  if (subscribe_mag && (topicMode == MagCalibrated)) {
    //  utilize magnetometer parameters
    for (int i = 0; i < 3; i++) {
      mm[i] -= magBias[i];
      mm[i] /= magScale[i];
    }
    eskf.setMagneticReference(magReference);
    eskf.setUsesMagnetometer(true);
  } else {
    //  fall back to simply gravity correction
    eskf.setUsesMagnetometer(false);
  }

  if (prevStamp.sec != 0) {
    double delta = imu->header.stamp.toSec() - prevStamp.toSec();
    eskf.predict(wm, delta);
  }
  prevStamp = imu->header.stamp;

  eskf.update(am, mm);

  const kr::quat<double> Q = eskf.getQuat(); //  updated quaternion
  const kr::AttitudeESKF::vec3 w =
      eskf.getAngularVelocity(); //  bias subtracted

  if (topicMode == MagInCalibration) {
    calib.appendSample(Q, mm);
    if (calib.isReady()) {

      try {
        calib.calibrate(kr::AttitudeMagCalib::FullCalibration);
      }
      catch (kr::AttitudeMagCalib::singular_hessian &e) {
        ROS_ERROR("Error: Failed to calibrate (error in minimization) : %s",
                  e.what());
        calib.reset();
      }

      if (calib.isCalibrated()) {
        auto b = calib.getBias();
        auto s = calib.getScale();
        auto r = calib.getReference();

        ROS_INFO("Calibration complete. Saving to rosparam.");
        if (preferredCalib == kr::AttitudeMagCalib::FullCalibration) {
          ROS_WARN("bias: %f, %f, %f", b[0], b[1], b[2]);
          ROS_WARN("scale: %f, %f, %f", s[0], s[1], s[2]);
          magScale = s;
          magBias = b;
        }

        ROS_WARN("ref: %f, %f, %f", r[0], r[1], r[2]);
        magReference = r;

        //  write to rosparam server
        for (int i = 0; i < 3; i++) {
          if (preferredCalib == kr::AttitudeMagCalib::FullCalibration) {
            // only write if full calibration requested
            nh->setParam("mag_calib/bias/" + sfx[i], magBias[i]);
            nh->setParam("mag_calib/scale/" + sfx[i], magScale[i]);
          }
          nh->setParam("mag_calib/reference/" + sfx[i], magReference[i]);
        }

        topicMode = MagCalibrated;
      }
    }
  }

  //  publish IMU topic
  sensor_msgs::Imu filtImu;
  filtImu.header.stamp = imu->header.stamp;
  filtImu.header.frame_id = bodyFrameName;

  filtImu.linear_acceleration = imu->linear_acceleration;
  filtImu.linear_acceleration_covariance = imu->linear_acceleration_covariance;
  filtImu.angular_velocity_covariance = imu->angular_velocity_covariance;

  tf::vectorEigenToMsg(w, filtImu.angular_velocity);

  filtImu.orientation.w = Q.w();
  filtImu.orientation.x = Q.x();
  filtImu.orientation.y = Q.y();
  filtImu.orientation.z = Q.z();

  //  append our covariance estimate to the new IMU message
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      filtImu.orientation_covariance[i * 3 + j] = eskf.getCovariance()(i, j);
    }
  }
  pubImu.publish(filtImu);

  //  publish bias
  geometry_msgs::Vector3Stamped bias;
  bias.header.stamp = filtImu.header.stamp;
  bias.header.frame_id = bodyFrameName;
  tf::vectorEigenToMsg(eskf.getGyroBias(), bias.vector);
  pubBias.publish(bias);

  //  publish status
  kr_attitude_eskf::Status status;
  status.header.stamp = filtImu.header.stamp;
  status.header.frame_id = "0";
  status.magStatus = topicMode;
  pubStatus.publish(status);

  //  publish compensated magnetic field
  if (subscribe_mag) {
    sensor_msgs::MagneticField fieldOut;
    fieldOut.header.stamp = filtImu.header.stamp;
    fieldOut.header.frame_id = bodyFrameName;
    tf::vectorEigenToMsg(mm, fieldOut.magnetic_field);
    pubField.publish(fieldOut);
  }

  if (publish_pose) {
    geometry_msgs::PoseStamped pose;
    pose.header.stamp = filtImu.header.stamp;
    pose.pose.orientation = filtImu.orientation;
    pose.pose.position.x = pose.pose.position.y = pose.pose.position.z = 0;
    pose.header.frame_id = bodyFrameName;
    pubPose.publish(pose);
  }
  
  //  broadcast frame
  if (broadcast_frame) {
    tf::Transform transform;
    transform.setRotation(tf::Quaternion(Q.x(), Q.y(), Q.z(), Q.w()));
    transform.setOrigin(tf::Vector3(0.0, 0.0, 0.0));

    tfBroadcaster->sendTransform(tf::StampedTransform(
        transform, filtImu.header.stamp, "/world", bodyFrameName));
  }
}

int main(int argc, char **argv) {
  ros::init(argc, argv, "kr_attitude_eskf");
  nh = ros::NodeHandlePtr(new ros::NodeHandle("~"));

  //  synchronized subscribers
  message_filters::Subscriber<sensor_msgs::Imu> imuSub;
  message_filters::Subscriber<sensor_msgs::MagneticField> fieldSub;
  message_filters::TimeSynchronizer<
      sensor_msgs::Imu, sensor_msgs::MagneticField> sync(imuSub, fieldSub, 1);
  ros::Subscriber imuSingleSub;

  //  find which topics to subscribe to
  std::string imu_topic="imu", field_topic="field";
  std::string calibrate_mag;
  bool calib_requested = false;

  nh->param("enable_magnetometer", subscribe_mag, false);
  nh->param("calibrate_magnetometer", calibrate_mag, std::string("none"));

  std::transform(calibrate_mag.begin(), calibrate_mag.end(),
                 calibrate_mag.begin(), ::tolower);
  if (calibrate_mag == "none") {
    calib_requested = false;
  } else if (calibrate_mag == "reference") {
    preferredCalib = kr::AttitudeMagCalib::ReferenceCalibration;
    calib_requested = true;
  } else if (calibrate_mag == "full") {
    preferredCalib = kr::AttitudeMagCalib::FullCalibration;
    calib_requested = true;
  } else {
    ROS_ERROR("Invalid parameter specified for calibrate_magnetometer: %s",
              calibrate_mag.c_str());
    return -1; // exit
  }

  if (calibrate_mag != "none" && !subscribe_mag) {
    ROS_ERROR("Invalid parameter pair: enable_magnetometer must be true if "
              "calibrate_magnetometer is true!");
    return -1; // exit
  }

  ROS_INFO("Subscribing to IMU: ~%s", imu_topic.c_str());

  //  filtered IMU output
  pubImu = nh->advertise<sensor_msgs::Imu>("filtered_imu", 1);
  pubBias = nh->advertise<geometry_msgs::Vector3Stamped>("bias", 1);
  pubStatus = nh->advertise<kr_attitude_eskf::Status>("status", 1);

  try {
    if (subscribe_mag) {
      ROS_INFO("Subscribing to magnetic field: ~%s", field_topic.c_str());

      //  subscribe to indicated topics using tight timing
      imuSub.subscribe(*nh, imu_topic, 20);
      fieldSub.subscribe(*nh, field_topic, 20);
      sync.registerCallback(boost::bind(&imu_callback, _1, _2));

      pubField = nh->advertise<sensor_msgs::MagneticField>("adjusted_field", 1);
    } else {
      //  only subscribe to IMU
      imuSingleSub = nh->subscribe<sensor_msgs::Imu>(
          imu_topic, 5,
          boost::bind(&imu_callback, _1, sensor_msgs::MagneticFieldConstPtr()));
    }
  }
  catch (ros::InvalidNameException &e) {
    //  invalid topic selected
    ROS_ERROR("ROS Exception (invalid topic): %s", e.what());
    return -1;
  }

  kr::AttitudeESKF::VarSettings var;
  double gyro_bias_thresh;

  //  load all remaining parameters
  nh->param("broadcast_frame", broadcast_frame, false);
  nh->param("publish_pose", publish_pose, false);
  
  if (publish_pose) {
    pubPose = nh->advertise<geometry_msgs::PoseStamped>("pose", 1);
  }
  
  //  noise parameters
  for (int i = 0; i < 3; i++) {
    nh->param("noise_std/accel/" + sfx[i], var.accel[i], 1.0);
    nh->param("noise_std/gyro/" + sfx[i], var.gyro[i], 0.001);
    nh->param("noise_std/mag/" + sfx[i], var.mag[i], 0.1);
  }
  nh->param("gyro_bias_thresh", gyro_bias_thresh, 1e-2); //  rad/s

  if (subscribe_mag) {
    if (calib_requested) {
      topicMode = MagInCalibration;
      ROS_WARN("Attitude ESKF loaded in calibration mode. Turn the device "
               "about all three axes!");
    } else {
      if (nh->hasParam("mag_calib/reference")) {
        ROS_INFO("Will load magnetometer reference from rosparam");
        topicMode = MagCalibrated;
      } else {
        ROS_WARN("Warning: No reference vector found. Did you calibrate?");
        topicMode = MagIdle;
      }
    }

    //  load magnetometer calibration, if set
    for (int i = 0; i < 3; i++) {
      nh->param("mag_calib/bias/" + sfx[i], magBias[i], 0.0);
      nh->param("mag_calib/scale/" + sfx[i], magScale[i], 1.0);
      nh->param("mag_calib/reference/" + sfx[i], magReference[i], 0.0);
    }
  }

  //  configure noise parameters
  eskf.setVariances(var);
  eskf.setEstimatesBias(true);
  eskf.setGyroBiasThreshold(gyro_bias_thresh);

  //  frame_id used for measurements
  bodyFrameName = ros::this_node::getName() + "/bodyFrame";
  
  if (broadcast_frame) {
    tfBroadcaster = boost::shared_ptr<tf::TransformBroadcaster>(
        new tf::TransformBroadcaster());
  }
 
  ros::spin();
  return 0;
}
