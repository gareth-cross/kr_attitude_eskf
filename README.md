# kr_attitude_eskf

kr_attitude_eskf is an implementation of the Error-State Kalman Filter described in:

* “Attitude Error Representations for Kalman Filtering” - F. Landis Markley

## Version History

* **0.0.8**:
 - Removed the reference calib option, which produced somewhat poor results anyways.
 - Mag calib code now uses either float or double.
 - Removed kr_math, eigen_conversions dependencies.
* **0.0.7**:
  - AttitudeESKF accepts matrices for covariance parameters, instead of assuming diagonal noise.
  - Noise parameters are not provided by `rosparam` anymore, and are expected to arrive in the incoming messages.
  - Process model treats gyro covariance with the correct units (rad/s)^2 instead of (rad)^2.
  - Magnetic bias/scale/reference vectors are loaded as lists instead of key-value arrays.
  - Magnetometer calibration now uses least squares for initial guess, and LM iteration is robustified with Cauchy weighting.
  - Calibration is triggered by a callable ROS service, instead of from the launch file options.
  - Added support for `diagnostic_updater`.
  - Renamed `adjusted_field` topic to `corrected_field`.
  - Refactored ROS code into a Node class.
  - Removed `broadcast_frame` option and `tf`-dependent code.
  - Pose is always published.
* **0.0.6**:
  - Added `publish_pose` option.
* **0.0.5**:
  - First public release.

## Description

This implementation utilizes accelerometer, gyroscope, and (optionally) the magnetometer in order to provide an attitude estimate. Gyroscope biases are estimated online with a moving average filter.

A 4th order Runge Kutta process model is employed, along with a classical EKF update - save for the error-state modifications. Internally, the filter uses a quaternion parameterization of SO(3), thereby avoiding the uglier qualities of euler angles.

## Dependencies

The core ESKF code depends on Eigen. The ROS node + magnetometer calibrator depends on [kr_math](https://github.com/KumarRobotics/kr_math) also.

## Parameters

The ROS implementation exposes several parameters:

|Parameter|Definition|Default|
|---|---|---|
|`enable_magnetometer`|If true, magnetometer readings are included in the state estimation update.|`false`|
|`mag_calib/bias`|Bias of the magnetometer.|`[0,0,0]`|
|`mag_calib/scale`|Scale of the magnetometer.|`[1,1,1]`|
|`mag_calib/reference`|World frame reference vector of the magnetometer.|`[0,0,0]`|
|`process_scale_factor`|'Fudge factor' to multiply by incoming gyro covariances.|1|
|`gyro_bias_thresh`|Threshold of motion below which we may estimate gyro bias.|0.01 rad/s|

When using the node, you should remap `~imu` and `~field` to the appropriate topics. See `attitude_eskf.launch` for an example.

## Topics

The following topics are always published:

|Name|Description|
|---|---|
|`filtered_imu`|Filter orientation, covariance, IMU acceleration, and angular velocity (with bias subtracted).|
|`bias`|Current estimate of gyroscope bias.|
|`status`|Status of filter. See `Status.msg` for details.|
|`pose`|Stamped pose, with translation component zeroed.|

If `enable_magnetometer` is set, the following topics will also be published:

|Name|Description|
|---|---|
|`corrected_field`|Magnetic field, with bias and scale estimates applied.|

## Calibration

This node includes a built-in magnetometer calibration mode. Calibration is started by calling the `begin_calibration` service. You will observe a log message on `rosout` indicating that calibration has started. Steps:

1. Rotate the device about the world z axis 360 degrees.
2. Rotate the device 90 degrees about the roll/pitch axis, and then again 360 degrees about the world z axis.
3. You should see the resulting parameters printed on `rosout`. They will also be saved to `rosparam`.
4. After calibration, `enable_magnetometer` will be set to true and the magnetic field will be included in the estimate.

## Important Operating Notes

* It is assumed that all IMU topics are **strictly synchronized**.
* It is assumed that all IMU sensors are **physically aligned**.
* The default value of `mag_calib/reference` is **not** valid. You must either calibrate the device, or set a known value using a magnetic model (such as [IGRF11](http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)).

## Bug?

Please report issues to the official package maintainer.

## Links

* [Early version running on iOS](http://www.youtube.com/watch?v=ijK2ndEGBXA)
* Avik De has used this filter on STM32 w/ MPU6000: [Details](http://avikde.me/IMU-Filtering-STM32-MPU6000/)
* [Video of a pico-quadrotor using ESKF](http://www.youtube.com/watch?v=HpO8DpeYC2k&spfreload=10)

