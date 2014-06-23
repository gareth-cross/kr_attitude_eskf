# kr_attitude_eskf

kr_attitude_eskf is an implementation of the Error-State Kalman Filter described
in:

* “Attitude Error Representations for Kalman Filtering” - F. Landis Markley

## Description

This implementation utilizes accelerometer, gyroscope, and (optionally) the magnetometer in order to provide an attitude estimate. Gyroscope biases are estimated online.

A 4th order Runge Kutta process model is employed, along with a classical EKF update - save for the error-state modifications. Internally, the filter uses a quaternion parameterization of SO(3), thereby avoiding the ugly qualities of euler angles.

## Dependencies

The core ESKF code depends only on `kr_math`. The ROS node depends on various other packages, listed in the manifest.

## Links

[Early version running on iOS](http://www.youtube.com/watch?v=ijK2ndEGBXA)
