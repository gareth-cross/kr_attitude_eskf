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

#include <eigen_conversions/eigen_msg.h>

#include <message_filters/subscriber.h>
#include <message_filters/time_synchronizer.h>

#include <tf/transform_broadcaster.h>
#include <tf/transform_datatypes.h>

#include <kr_attitude_eskf/BeginCalibration.h>
#include <kr_attitude_eskf/Status.h>

#include "AttitudeESKF.hpp"
#include "SO3.hpp"

using namespace std;
using namespace Eigen;

kr::AttitudeESKF eskf;
ros::Publisher pubImu;
ros::Publisher pubBias;
ros::Publisher pubField;
ros::Publisher pubStatus;

ros::NodeHandle* cur_nh;  //  temp hack

const std::string sfx[] = {"x", "y", "z"};  //  suffixes

boost::shared_ptr<tf::TransformBroadcaster> tfBroadcaster;  //  for broadcasting body frame
std::string bodyFrameName;
bool broadcast_frame = false;

enum {
  MagIdle=0,
  MagInCalibration=1,
  MagCalibrated=2,
} topicMode = MagIdle;
bool subscribe_mag = false; //  use the magnetometer or not?

#define BIN_DIM  (40)

struct SampleBin {
  Vector3d field; //  field measured in this sample
  kr::quat<double> q; //  unreferenced quaternion for this sample
};

std::map<int, SampleBin> binsYaw;
std::map<int, SampleBin> binsPitch;
std::map<int, SampleBin> binsRoll;

//  default bias and scale
Vector3d magBias = Vector3d::Zero();
Vector3d magScale = Vector3d::Ones();
Vector3d magReference = Vector3d::Zero();

bool begin_calibration(kr_attitude_eskf::BeginCalibration::Request& req,
                       kr_attitude_eskf::BeginCalibration::Response& resp)
{
  if (topicMode == MagInCalibration) {
    resp.errorCode = -1;  //  already calibrating
  } else {
    topicMode = MagInCalibration;
    resp.errorCode = 0; //  transition to calibration mode
  }

  return true;
}

void imu_callback(const sensor_msgs::ImuConstPtr& imu, const sensor_msgs::MagneticFieldConstPtr& field)
{
  Vector3d wm;  //  measured angular rate
  Vector3d am;  //  measured acceleration
  Vector3d mm;  //  measured magnetic field

  tf::vectorMsgToEigen(imu->angular_velocity, wm);
  tf::vectorMsgToEigen(imu->linear_acceleration, am);

  if (subscribe_mag) {
    tf::vectorMsgToEigen(field->magnetic_field, mm);
  } else {
    mm.setZero(); //  safe default
  }

  if (subscribe_mag && (topicMode == MagCalibrated))
  {
    //  utilize magnetometer parameters
    for (int i=0; i < 3; i++) {
      mm[i] -= magBias[i];
      mm[i] /= magScale[i];
    }
    eskf.setMagneticReference(magReference);
    eskf.setUsesMagnetometer(true);
  }
  else {
    //  fall back to simply gravity correction
    eskf.setUsesMagnetometer(false);
  }

  eskf.predict(wm,imu->header.stamp.toSec());
  eskf.update(am,mm);

  kr::quat<double> Q = eskf.getQuat();                    //  updated quaternion
  kr::AttitudeESKF::vec3 w = eskf.getAngularVelocity();   //  bias subtracted

  if (topicMode == MagInCalibration)
  {
    Vector3d angs = kr::getRPY(Q.to_matrix());

    int idx_roll = std::floor( (angs[0] + M_PI) / (2*M_PI) * BIN_DIM );
    int idx_pitch = std::floor( (angs[1] + M_PI/2) / M_PI * BIN_DIM );
    int idx_yaw = std::floor( (angs[2] + M_PI) / (2*M_PI) * BIN_DIM );

    {
      SampleBin bin;
      bin.field = mm;
      bin.q = Q;

      binsRoll[idx_roll] = bin;
      binsYaw[idx_yaw] = bin;
      binsPitch[idx_pitch] = bin;

      ROS_INFO("R: %lu, P: %lu, Y: %lu", binsRoll.size(), binsPitch.size(), binsYaw.size());

      bool enough = (binsRoll.size() > BIN_DIM*0.9) && (binsPitch.size() > BIN_DIM*0.9) && (binsYaw.size() > BIN_DIM*0.9);

      if ( enough )
      {
        std::vector<SampleBin> sampleBins;
        for (auto it : binsRoll) {
          sampleBins.push_back(it.second);
        }
        for (auto it : binsYaw) {
          sampleBins.push_back(it.second);
        }
        for (auto it : binsPitch) {
          sampleBins.push_back(it.second);
        }

        ROS_INFO("Collected enough bins to calibrate");

        vector<Vector3d> samples;

        Vector3d max,min;
        min[0] = min[1] = min[2] = std::numeric_limits<double>::infinity();
        max[0] = max[1] = max[2] = -min[0]; //  -infinity

        double mean_rad = 0.0;
        double mean_rad_sqr = 0.0;

        for (auto i = sampleBins.begin(); i != sampleBins.end(); i++)
        {
          samples.push_back(i->field);

          const double r = i->field.norm();
          mean_rad += r;
          mean_rad_sqr += r*r;
        }
        mean_rad /= sampleBins.size();      //  mean radius
        mean_rad_sqr /= sampleBins.size();

        //  standard deviation
        const float std_dev = std::sqrt(mean_rad_sqr - mean_rad*mean_rad);

        //  calculate bias term
        for (Vector3d& s : samples) {
          for (int j=0; j < 3; j++) {
            max[j] = std::max(max[j], s[j]);
            min[j] = std::min(min[j], s[j]);
          }
        }

        //  estimate of soft bias
        Vector3d bias = (max + min) / 2;

        //  reject outliers
        auto j = sampleBins.begin();
        for (auto i = samples.begin(); i != samples.end();)
        {
          const Vector3d s = *i - bias;

          if (std::abs(s.norm() - mean_rad) > std_dev*2) {
            i = samples.erase(i);
            j = sampleBins.erase(j);
          } else {
            i++;
            j++;
          }
        }

        ROS_INFO("%lu samples left after discarding outliers", samples.size());
        ROS_INFO("bias: %f, %f, %f", bias[0], bias[1], bias[2]);

        //  attempt to determine scale factors via GN-NLS
        Vector3d scl;
        scl.setOnes();

        int num_iter = 0;
        while (num_iter++ < 10)
        {
          MatrixXd J(samples.size(), 6);
          VectorXd r(samples.size());

          for (size_t i=0; i < samples.size(); i++)
          {
            double x = (samples[i][0] - bias[0]) / scl[0];
            double y = (samples[i][1] - bias[1]) / scl[1];
            double z = (samples[i][2] - bias[2]) / scl[2];

            double rad2 = x*x + y*y + z*z;
            r[i] = mean_rad*mean_rad - rad2;

            J(i,0) = -2 * x*x / scl[0];
            J(i,1) = -2 * y*y / scl[1];
            J(i,2) = -2 * z*z / scl[2];

            J(i,3) = -2 * x / scl[0];
            J(i,4) = -2 * y / scl[1];
            J(i,5) = -2 * z / scl[2];
          }

          Matrix<double,6,6> H = J.transpose() * J;
          for (int i=0; i < H.rows(); i++) {
            H(i,i) *= 1.01;
          }

          bool invertible;
          decltype(H) Hinv;

          auto LU = H.fullPivLu();
          invertible = LU.isInvertible();
          Hinv = LU.inverse();

          Matrix<double,6,1> update = Hinv * J.transpose() * r;

          bias += update.block<3,1>(3,0);
          scl += update.block<3,1>(0,0);

          if (!invertible) {
            ROS_ERROR("Failed to optimize");
            return;
          }
        }

        //bias.setZero();
        //scl.setOnes();

        magBias = bias;
        magScale = scl;

        ROS_INFO("Adjusted scale: %f, %f, %f", scl[0], scl[1], scl[2]);
        ROS_INFO("Adjusted bias: %f, %f, %f", bias[0], bias[1], bias[2]);

        double hAvg = 0.0;
        double vAvg = 0.0;
        long count=0;

        //  determine the inertial field
        for (auto i = sampleBins.begin(); i != sampleBins.end(); i++)
        {
          const SampleBin& bin = *i;

          Matrix3d wRb = bin.q.to_matrix().cast<double>();

          Vector3d acc = wRb.transpose().block<3,1>(0,2); //  body frame gravity vector
          if (acc[2] < -0.8f)
          {
            //  close to the vertical, use this sample
            Vector3d rpy = kr::getRPY(wRb);

            //  rotate back to level
            Matrix3d trans = kr::rotation_y(rpy[1]) * kr::rotation_x(rpy[0]);
            Vector3d f = bin.field;

            for (int j=0; j < 3; j++) {
              f[j] -= bias[j];
              f[j] /= scl[j];
            }

            Vector3d mw = trans * f;

            double horz = std::sqrt(mw[0]*mw[0] + mw[1]*mw[1]);
            double vert = mw[2];

            hAvg += horz;
            vAvg += vert;
            count++;

            ROS_INFO("H: %f, V: %f", horz, vert);
          }
        }

        hAvg /= count;
        vAvg /= count;

        Eigen::Matrix<double,3,1> vRef;
        vRef.setZero();
        vRef[0] = hAvg;
        vRef[2] = vAvg;

        magReference = vRef;

        //  save to rosparam
        for (int i=0; i < 3; i++) {
          cur_nh->setParam("mag_calib/bias/"+sfx[i],magBias[i]);
          cur_nh->setParam("mag_calib/scale/"+sfx[i],magScale[i]);
          cur_nh->setParam("mag_calib/reference/"+sfx[i],magReference[i]);
        }

        eskf.setMagneticReference(vRef);
        topicMode = MagCalibrated;
      }
    }
  }

  //  publish IMU topic
  sensor_msgs::Imu filtImu;
  filtImu.header.stamp = ros::Time::now();

  filtImu.linear_acceleration = imu->linear_acceleration;
  filtImu.linear_acceleration_covariance = imu->linear_acceleration_covariance;
  filtImu.angular_velocity_covariance = imu->angular_velocity_covariance;

  filtImu.angular_velocity.x = w[0];
  filtImu.angular_velocity.y = w[1];
  filtImu.angular_velocity.z = w[2];

  filtImu.orientation.w = Q.a();
  filtImu.orientation.x = Q.b();
  filtImu.orientation.y = Q.c();
  filtImu.orientation.z = Q.d();

  //  append our covariance estimate to the new IMU message
  for (int i=0; i < 3; i++) {
    for (int j=0; j < 3; j++) {
      filtImu.orientation_covariance[i*3 + j] = eskf.getCovariance()(i,j);
    }
  }
  pubImu.publish(filtImu);

  //  publish bias
  geometry_msgs::Vector3Stamped bias;
  bias.header.stamp = filtImu.header.stamp;
  bias.vector.x = eskf.getGyroBias()[0];
  bias.vector.y = eskf.getGyroBias()[1];
  bias.vector.z = eskf.getGyroBias()[2];
  pubBias.publish(bias);

  //  publish status
  kr_attitude_eskf::Status status;
  status.header.stamp = filtImu.header.stamp;
  status.magStatus = topicMode;
  pubStatus.publish(status);

  //  publish compensated magnetic field
  if (subscribe_mag) {
    sensor_msgs::MagneticField fieldOut;
    fieldOut.header.stamp = filtImu.header.stamp;
    fieldOut.magnetic_field.x = mm[0];
    fieldOut.magnetic_field.y = mm[1];
    fieldOut.magnetic_field.z = mm[2];
    pubField.publish(fieldOut);
  }

  //  broadcast frame
  if (broadcast_frame) {
    tf::Transform transform;
    transform.setRotation( tf::Quaternion(Q.b(), Q.c(), Q.d(), Q.a()) );
    transform.setOrigin( tf::Vector3(0.0, 0.0, 0.0) );

    tfBroadcaster->sendTransform( tf::StampedTransform(transform, filtImu.header.stamp, "fixedFrame", bodyFrameName) );
  }
}

int main(int argc, char ** argv)
{
  ros::init(argc, argv, "kr_attitude_eskf");
  ros::NodeHandle nh("~");
  ros::NodeHandle nh_pub;  //  public
  ros::ServiceServer calibSrv;
  cur_nh = &nh;

  //  synchronized subscribers
  message_filters::Subscriber<sensor_msgs::Imu> imuSub;
  message_filters::Subscriber<sensor_msgs::MagneticField> fieldSub;
  message_filters::TimeSynchronizer<sensor_msgs::Imu, sensor_msgs::MagneticField> sync(imuSub,fieldSub,1);

  //  find which topics to subscribe to
  std::string imu_topic, field_topic;
  nh.param("imu_topic", imu_topic, std::string("/imu/imu"));
  nh.param("field_topic", field_topic, std::string("/imu/magnetic_field"));
  nh.param("enable_magnetometer", subscribe_mag, false);

  ROS_INFO("Subscribing to IMU: %s", imu_topic.c_str());

  //  filtered IMU output
  pubImu = nh.advertise<sensor_msgs::Imu>("filtered_imu", 1);
  pubBias = nh.advertise<geometry_msgs::Vector3Stamped>("bias", 1);
  pubStatus = nh.advertise<kr_attitude_eskf::Status>("status", 1);

  try
  {
    if (subscribe_mag)
    {
      ROS_INFO("Subscribing to magnetic field: %s", field_topic.c_str());

      //  subscribe to indicated topics using tight timing
      imuSub.subscribe(nh,imu_topic,20);
      fieldSub.subscribe(nh,field_topic,20);
      sync.registerCallback(boost::bind(&imu_callback, _1, _2));

      pubField = nh.advertise<sensor_msgs::MagneticField>("adjusted_field", 1);

      //  service for triggering calibration
      calibSrv = nh_pub.advertiseService(ros::this_node::getName() + "/begin_calibration", begin_calibration);
    }
    else
    {
      //  only subscribe to IMU
      nh.subscribe<sensor_msgs::Imu>(imu_topic, 5, boost::bind(&imu_callback, _1, sensor_msgs::MagneticFieldConstPtr()));
    }
  }
  catch(ros::InvalidNameException& e) {
    //  invalid topic selected
    ROS_ERROR("ROS Exception (invalid topic): %s", e.what());
    return -1;
  }

  kr::AttitudeESKF::VarSettings var;
  double gyro_bias_thresh;

  //  load all parameters
  nh.param("broadcast_frame", broadcast_frame, false);

  //  noise parameters
  for (int i=0; i < 3; i++) {
    nh.param("noise_std/accel/"+sfx[i], var.accel[i], 1.0);
    nh.param("noise_std/gyro/"+sfx[i], var.gyro[i], 0.001);
    nh.param("noise_std/mag/"+sfx[i], var.mag[i], 0.1);
  }
  nh.param("gyro_bias_thresh", gyro_bias_thresh, 1e-2); //  rad/s

  if (subscribe_mag)
  {
    if (nh.hasParam("mag_calib/bias") &&
        nh.hasParam("mag_calib/scale") &&
        nh.hasParam("mag_calib/reference"))
    {
      ROS_INFO("Loading mag bias and scale terms from rosparam");
      topicMode = MagCalibrated;
    }
    else {
      ROS_WARN("Warning: No calibration parameters found. Did you load your yaml file?");
      topicMode = MagIdle;
    }

    //  magnetometer calibration
    for (int i=0; i < 3; i++) {
      nh.param("mag_calib/bias/"+sfx[i], magBias[i], 0.0);
      nh.param("mag_calib/scale/"+sfx[i], magScale[i], 1.0);
      nh.param("mag_calib/reference/"+sfx[i], magReference[i], 0.0);
    }
  }

  eskf.setVariances(var);
  eskf.setEstimatesBias(true);
  eskf.setGyroBiasThreshold(gyro_bias_thresh);

  if (broadcast_frame)
  {
    tfBroadcaster = boost::shared_ptr<tf::TransformBroadcaster>( new tf::TransformBroadcaster() );
    bodyFrameName = ros::this_node::getName() + "/bodyFrame";
  }

  ros::spin();
  return 0;
}
