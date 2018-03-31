# Linear_Kalman_Filter

Adapted from Simon D. Levy's [TinyEKF](https://github.com/simondlevy/TinyEKF). I mainly changed the tiny_ekf.c file to make the filter works with linear models. Current application is for yaw angle estimation, still working on testing with Teensy Arduino; still working on documentation. The interface should be able to work with different models. 

To understand how Kalman filter works, see [Bzarg's blogpost with illustrations](http://www.bzarg.com/p/how-a-kalman-filter-works-in-pictures/) and [Simon's interactive tutorial](https://home.wlu.edu/~levys/kalman_tutorial/)