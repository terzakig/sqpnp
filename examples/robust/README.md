
This directory demonstrates [SQPnP](https://github.com/terzakig/sqpnp) with the [RansacLib](https://github.com/tsattler/RansacLib) library.

## General
The robust_sqpnp_example demo program reads from two files: one with the camera intrinsics K and another containing the 3D - 2D corresponding points in separate lines as<br>
X Y Z x y<br>
where x y are the image projections in pixels. To use normalized coordinates for x y, specify either the identity matrix in K or - as its filename.

For more advanced pose estimation (e.g., a non-linear refinement on the inliers, unknown focal length, etc), check the [posest](https://users.ics.forth.gr/~lourakis/posest/) library.

## Running
./robust_sqpnp_example  data/K.txt  data/32D.txt

## Caveat
RANSAC parameters such as the inlier percentage, outlier threshold, number of iterations, etc should be tuned to your particular problem!
