# kr_attitude_eskf

kr_attitude_eskf is an implementation of the Error-State Kalman Filter described
in:

* “Attitude Error Representations for Kalman Filtering” - F. Landis Markley

## Description

This implementation utilizes accelerometer, gyroscope, and (optionally) the magnetometer in order to provide an attitude estimate. Gyroscope biases are estimated online.

A 4th order Runge Kutta process model is employed, along with a classical EKF update - save for the error-state modifications. Internally, the filter uses a quaternion parameterization of SO(3), thereby avoiding the uglier qualities of euler angles.

## Dependencies

The core ESKF code depends on Eigen. The ROS node + magnetometer calibrator depends on [kr_math](https://github.com/KumarRobotics/kr_math) also.

## Parameters

The ROS implementation exposes several parameters:

|Parameter|Definition|Default|
|---|---|---|
|`enable_magnetometer`|If true, magnetometer readings are included in the state estimation update.|`false`|
|`calibrate_magnetometer`|Type of calibration - see below.|`"none"`|
|`broadcast_frame`|If true, the TF frame for the body to world transform is broadcast.|`false`|
|`noise_std/accel/[x,y,z]`|Variance of the accelerometer noise.|`0.1`|
|`noise_std/gyro/[x,y,z]`|Variance of the gyro noise.|`0.01`|
|`noise_std/mag/[x,y,z]`|Variance of the magnetometer noise.|`0.1`|
|`mag_calib/bias/[x,y,z]`|Bias of the magnetometer.|`0`|
|`mag_calib/scale/[x,y,z]`|Scale of the magnetometer.|`1`|
|`mag_calib/reference/[x,y,z]`|World frame reference vector of the magnetometer.|`0`|

When using the node, you should remap `~imu` and `~field` to the appropriate topics. See `attitude_eskf.launch` for an example.

## Topics

The following topics are always published:

|Name|Description|
|---|---|
|`filtered_imu`|Filter orientation, covariance, IMU acceleration, and angular velocity (with bias subtracted).|
|`bias`|Current estimate of gyroscope bias.|
|`status`|Status of filter. See `Status.msg` for details.|

If `enable_magnetometer` is set, the following topics will also be published:

|Name|Description|
|---|---|
|`adjusted_field`|Magnetic field, with bias and scale estimates applied.|

## Calibration

This node includes a built-in magnetometer calibration mode. Calibration modes are selected with the `calibrate_magnetometer` option, which accepts the following strings:

* `none`: No calibration is performed (default).
* `reference`: Only the magnetic reference vector (north) is determined.
* `full`: Bias, scale and reference vectors are determined.

When launching in calibration mode, it is expected that you will rotate the device about all three of its major axes. The node will display the resulting parameters when calibration is complete. **The calculated parameters will also be saved to rosparam.**

Calibration is performed using non-linear weighted least squares.

## Important Operating Notes

* It is assumed that **all IMU topics are strictly synchronized**.
* It is assumed that **all IMU sensors are physically aligned**.
* If `calibrate_magnetometer` is set, then `enable_magnetometer` must also be set.
* The default value of `mag_calib/reference` is **not** valid. You must either calibrate the device, or set a known value using a magnetic model (such as [IGRF11](http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)).

## Bug?

Please report issues to the official package maintainer.

## Links

Demo: [Early version running on iOS](http://www.youtube.com/watch?v=ijK2ndEGBXA)
