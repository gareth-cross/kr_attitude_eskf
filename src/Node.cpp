/*
 * Node.cpp
 *
 *  Copyright (c) 2014 Gareth Cross. All rights reserved.
 *
 *  This file is part of kr_attitude_eskf.
 *
 *	Created on: 31/08/2014
 *		  Author: gareth
 */

#include <kr_attitude_eskf/Node.hpp>
#include <geometry_msgs/Vector3Stamped.h>
#include <geometry_msgs/PoseStamped.h>
#include <eigen_conversions/eigen_msg.h>
#include <assert.h>

namespace kr_attitude_eskf {

Node::Node(const ros::NodeHandle &nh) : nh_(nh), 
  sync_(TimeSyncPolicy(kROSQueueSize), subImu_, subField_), 
  calibState_(Uncalibrated) {
 
  //  load settings
  nh_.param("enable_magnetometer", enableMag_, false);
  nh_.param("gyro_bias_thresh",  gyroBiasThresh_, 0.01);
  //  load magnetometer calibration
  magBias_.setZero();
  magScale_.setConstant(1);
  magReference_.setZero();
  if (enableMag_) {
    try {
      std::vector<double> biasVec, scaleVec, refVec;
      nh_.getParam("mag_calib/bias", biasVec);
      nh_.getParam("mag_calib/scale", scaleVec);
      nh_.getParam("mag_calib/reference", refVec);
      if (biasVec.size() != 3 || scaleVec.size() != 3 || refVec.size() != 3) {
        throw std::invalid_argument("List w/ invalid size in mag_calib.");
      }
      
      magBias_ << biasVec[0], biasVec[1], biasVec[2];
      magScale_ << scaleVec[0], scaleVec[1], scaleVec[2];
      magReference_ << refVec[0], refVec[1], refVec[2];
      
      ROS_INFO("Loaded mag_calib:");
      ROS_INFO("Bias: %f,%f,%f",magBias_[0],magBias_[1],magBias_[2]);
      ROS_INFO("Scale: %f,%f,%f",magScale_[0],magScale_[1],magScale_[2]);
      ROS_INFO("Reference: %f,%f,%f",magReference_[0],magReference_[1],
          magReference_[2]);
      
      calibState_ = CalibrationComplete;
    } catch (std::exception& e) {
      ROS_WARN("Warning: Failed to load mag_calib parameters.");
      ROS_WARN("Reason: %s", e.what());
    }
    
    //  subscribe using time sync
    ROS_INFO("Subscribing to ~imu and ~magnetic_field.");
    
    subImu_.subscribe(nh_,"imu",kROSQueueSize);
    subField_.subscribe(nh_,"magnetic_field",kROSQueueSize);
    sync_.registerCallback(boost::bind(&Node::inputCallback, this, _1, _2));
  } else {
    //  subscribe only to IMU
    ROS_INFO("Subscribing to ~imu.");
    subImuUnsync_ = nh_.subscribe<sensor_msgs::Imu>("imu", kROSQueueSize,
                         boost::bind(&Node::inputCallback, this, _1, 
                                     sensor_msgs::MagneticFieldConstPtr()));
  }
  
  if (enableMag_ && calibState_==Uncalibrated) {
    ROS_WARN("Warning: Magnetometer is enabled but not calibrated.");
  }
  
  pubImu_ = nh_.advertise<sensor_msgs::Imu>("filtered_imu", 1);
  pubBias_ = nh_.advertise<geometry_msgs::Vector3Stamped>("bias", 1);
  pubPose_ = nh_.advertise<geometry_msgs::PoseStamped>("pose", 1);
  
  //  configure filter
  eskf_.setEstimatesBias(true);
  eskf_.setGyroBiasThreshold(gyroBiasThresh_);
  eskf_.getCovariance().setIdentity();  //  ~60deg std. on each axis
}

void Node::inputCallback(const sensor_msgs::ImuConstPtr& imuMsg,
                         const sensor_msgs::MagneticFieldConstPtr& magMsg) {
  if (enableMag_) {
    assert(magMsg);
  }
  
  kr::vec3d wm; //  measured angular rate
  kr::vec3d am; //  measured acceleration
  kr::vec3d mm; //  measured magnetic field

  tf::vectorMsgToEigen(imuMsg->angular_velocity, wm);
  tf::vectorMsgToEigen(imuMsg->linear_acceleration, am);
  if (magMsg) {
    tf::vectorMsgToEigen(magMsg->magnetic_field, mm);
  } else {
    mm.setZero();
  }
  //  copy covariances from message to matrix
  kr::AttitudeESKF::mat3 wCov, aCov, mCov = kr::AttitudeESKF::mat3::Zero();
  for (int i=0; i < 3; i++) {
    for (int j=0; j < 3; j++) {
      wCov(i,j) = imuMsg->angular_velocity_covariance[i*3 + j];
      aCov(i,j) = imuMsg->linear_acceleration_covariance[i*3 + j];
      if (magMsg) {
        mCov(i,j) = magMsg->magnetic_field_covariance[i*3 + j];
      }
    }
  }
  
  if (enableMag_ && calibState_==CalibrationComplete) {
    //  correct magnetic field
    mm.noalias() -= magBias_;
    for (int i=0; i < 3; i++) {
      mm[i] /= magScale_[i];
    }
    eskf_.setMagneticReference(magReference_);
    eskf_.setUsesMagnetometer(true);
  } else {
    eskf_.setUsesMagnetometer(false);
  }
  
  static ros::Time firstTime = imuMsg->header.stamp;
  if (prevStamp_.sec != 0) {
    //  run kalman filter
    const double delta = imuMsg->header.stamp.toSec() - prevStamp_.toSec();
    eskf_.predict(wm, delta, wCov);
    eskf_.update(am, aCov, mm, mCov);
    
    const kr::quatd wQb = eskf_.getQuat();          // updated quaternion
    const kr::vec3d w = eskf_.getAngularVelocity(); // ang vel. minus bias
   
    if (calibState_ == Calibrating || true) {
      //  update the calibrator   
      if ((imuMsg->header.stamp - firstTime).toSec() > 5) {
        calib_.appendSample(wQb, mm);
      }
      
//      if (calib_.isReady()) {
        
        
//        //  calibrate bias, scale and reference vector
//        calib_.calibrate(kr::AttitudeMagCalib::FullCalibration);
//      }
    }
    
    sensor_msgs::Imu imu = *imuMsg;
    imu.header.seq = 0;
    
    tf::vectorEigenToMsg(w,imu.angular_velocity);
    imu.orientation.w = wQb.w();
    imu.orientation.x = wQb.x();
    imu.orientation.y = wQb.y();
    imu.orientation.z = wQb.z();
    //  append our covariance estimate to the new IMU message
    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        imu.orientation_covariance[i*3 + j] = eskf_.getCovariance()(i,j);
      }
    }
    geometry_msgs::Vector3Stamped bias;
    bias.header = imu.header;
    tf::vectorEigenToMsg(eskf_.getGyroBias(),bias.vector);
    //  pose with no position
    geometry_msgs::PoseStamped pose;
    pose.header = imu.header;
    pose.pose.orientation = imu.orientation;
    
    pubImu_.publish(imu);
    pubBias_.publish(bias);
    pubPose_.publish(pose);
  }
  prevStamp_ = imuMsg->header.stamp;
}

}
