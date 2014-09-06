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
#include <kr_math/SO3.hpp>
#include <assert.h>

namespace kr_attitude_eskf {

Node::Node(const ros::NodeHandle &nh, const ros::NodeHandle &pnh) : nh_(pnh), 
  sync_(TimeSyncPolicy(kROSQueueSize), subImu_, subField_), 
  calibState_(Uncalibrated), init_(false) {
 
  //  load settings
  nh_.param("enable_magnetometer", enableMag_, false);
  nh_.param("gyro_bias_thresh", gyroBiasThresh_, 0.01);
  nh_.param("process_scale_factor", processScaleFactor_, 1.0);
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
  //  'unbiased' and scaled magnetic field
  pubField_ = nh_.advertise<sensor_msgs::MagneticField>("corrected_field",1);
  
  //  register updater
  if (!nh_.hasParam("diagnostic_period")) {
    nh_.setParam("diagnostic_period", 0.2);
  }
  updater_.setHardwareID(ros::this_node::getName());
  updater_.add("diagnostic_info", 
               boost::bind(&Node::diagnosticCallback, this, _1));
  
  //  register for service callback
  ros::NodeHandle publicNh(nh);
  srvCalibrate_ = publicNh.advertiseService(ros::this_node::getName() + "/begin_calibration", 
                                            &Node::beginCalibration, this);
  
  //  configure filter
  eskf_.setEstimatesBias(true);
  eskf_.setGyroBiasThreshold(gyroBiasThresh_);
}

void Node::saveCalibration() {
  std::vector<double> vec(3);
  vec[0] = magBias_[0];
  vec[1] = magBias_[1];
  vec[2] = magBias_[2];
  nh_.setParam("mag_calib/bias", vec);
  vec[0] = magScale_[0];
  vec[1] = magScale_[1];
  vec[2] = magScale_[2];
  nh_.setParam("mag_calib/scale", vec);
  vec[0] = magReference_[0];
  vec[1] = magReference_[1];
  vec[2] = magReference_[2];
  nh_.setParam("mag_calib/reference", vec);
}

void Node::inputCallback(const sensor_msgs::ImuConstPtr& imuMsg,
                         const sensor_msgs::MagneticFieldConstPtr& magMsg) {
  if (enableMag_) {
    assert(magMsg);
  }
  
  kr::vec3d wm; //  measured angular rate
  kr::vec3d am; //  measured acceleration
  kr::vec3d mm; //  measured magnetic field

  //std::cout << mm.transpose() << std::endl;
  
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
  wCov *= processScaleFactor_;  //  apply fudge factor
  
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
  
  if (prevStamp_.sec != 0) {
    if (!init_) {
      const kr::vec3d aCovDiag(aCov(0,0),aCov(1,1),aCov(2,2));
      const kr::vec3d mCovDiag(mCov(0,0),mCov(1,1),mCov(2,2));
      eskf_.initialize(am,aCovDiag,mm,mCovDiag);
      eskf_.getCovariance() = kr::mat3d::Identity()*5;  //  large uncertainty
      init_ = true;
      ROS_INFO("Initialized ESKF");
    }
    
    //  run kalman filter
    const double delta = imuMsg->header.stamp.toSec() - prevStamp_.toSec();
    eskf_.predict(wm, delta, wCov);
    eskf_.update(am, aCov, mm, mCov);
    
    const kr::quatd wQb = eskf_.getQuat();          // updated quaternion
    const kr::vec3d w = eskf_.getAngularVelocity(); // ang vel. minus bias
   
    if (calibState_ == Calibrating) {
      //  update the calibrator
      if (!calib_.isCalibrated()) {
        calib_.appendSample(wQb, mm);
        if (calib_.isReady()) {
          ROS_INFO("Collected sufficient samples. Calibrating...");
          //  calibrate bias, scale and reference vector
          try {
            calib_.calibrate(kr::AttitudeMagCalib::FullCalibration);
            
            ROS_INFO_STREAM("Bias: " << calib_.getBias().transpose());
            ROS_INFO_STREAM("Scale: " << calib_.getScale().transpose());
            ROS_INFO_STREAM("Reference: " << calib_.getReference().transpose());
            
            magBias_ = calib_.getBias();
            magScale_ = calib_.getScale();
            magReference_ = calib_.getReference();
            
            //  save to rosparam
            nh_.setParam("enable_magnetometer", true);
            saveCalibration();
            
            //  done, we can use this new calibration immediately
            calibState_ = CalibrationComplete;
            enableMag_ = true;
          }
          catch (std::exception& e) {
            ROS_ERROR("Calibration failed: %s", e.what());    
          }
        }
      }
    }

//    const kr::vec3d derot = eskf_.getPredictedField();
//    static ros::Publisher pubDebug = nh_.advertise<geometry_msgs::Vector3Stamped>("predicted", 1);
//    geometry_msgs::Vector3Stamped vec;
//    tf::vectorEigenToMsg(derot, vec.vector);
//    vec.header.stamp = imuMsg->header.stamp;
//    pubDebug.publish(vec);
    
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
    
    if (enableMag_) {
      sensor_msgs::MagneticField field = *magMsg;
      field.header.seq = 0;
      tf::vectorEigenToMsg(mm, field.magnetic_field);
      pubField_.publish(field);
    }
    
    pubImu_.publish(imu);
    pubBias_.publish(bias);
    pubPose_.publish(pose);
  }
  prevStamp_ = imuMsg->header.stamp;
  updater_.update();
}

bool Node::beginCalibration(std_srvs::Empty::Request&,
                            std_srvs::Empty::Response&) {
  //  enter calibration mode
  calibState_ = Calibrating;
  calib_.reset();
  ROS_INFO("Entering calibration mode.");
  return true;
}

void 
Node::diagnosticCallback(diagnostic_updater::DiagnosticStatusWrapper &stat) {
  
  stat.summary(diagnostic_msgs::DiagnosticStatus::OK, "Running");
  
  stat.add("Magnetometer enabled:", enableMag_);
  stat.add("Gyro bias threshold:", gyroBiasThresh_);
  stat.add("Process scale factor:", processScaleFactor_);
  std::stringstream ss;
  //  current bias estimate
  ss << eskf_.getGyroBias().transpose();
  stat.add("Gyroscope bias:", ss.str());
  //  mag calib state
  std::string calib;
  if (calibState_ == Uncalibrated) {
    calib = "Uncalibrated";
  } else if (calibState_ == Calibrating) {
    calib = "Calibrating";
  } else {
    calib = "Complete";
  }
  stat.add("Calibration state:", calib);
}

}
