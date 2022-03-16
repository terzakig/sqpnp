# SQPnP 
C++ Implementation of the SQPnP algorithm.

The algorithm is the generic PnP solver described in the paper ["A Consistently Fast and Globally Optimal Solution to the Perspective-n-Point Problem" by G. Terzakis and M. Lourakis](http://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460.pdf). For more intuition, refer to the supplementary material [here](https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460-supp.pdf).

## Required libraries
SQPnP requires the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library to build. Besides rank revealing QR and optionally SVD, the use of Eigen is confined to matrix addition, transposition and multiplication.
Choosing Eigen was motivated by its increasing popularity and lightweight character. There are also three examples of using the solver in this repository:
1. ) one that uses OpenCV, just for the sake of demonstrating the initialization of SQPnP with ``cv::Point_<>`` and ``cv::Point3_<>`` structures,
2. ) another in which data points are copied from plain 2D arrays, and
3. ) a third which demonstrates SQPnP within RANSAC.

Build will proceed with either one of 1) or 2), depending on whether OpenCV is found or not.
Example 3) uses the [RansacLib](https://github.com/tsattler/RansacLib) template-based library, (part of) which is included in this repository.

## Build
-----

Create a ``build`` directory in the root of the cloned repository and run ``cmake``:

``mkdir build``

``cd build``

``cmake ..``

or, for a *release* build,

``cmake .. -DCMAKE_BUILD_TYPE=Release``

The latter will allow for more accurate timing of average execution time. Finally build everything:

``make``

To run the PnP example, once in the ``build`` directory,

``./examples/sqpnp_example``

To run the robust estimation example, from the ``build`` directory,
``./examples/robust_sqpnp_example data/K.txt data/32D.txt``

## Non-default parameters
See ``struct SolverParameters`` in ``types.h`` which contains SQPnP's parameters that can be specified by the caller.
For instance, to use SVD instead of RRQR for the nullspace basis of Omega, the following fragment can be used:
```c++
  // call solver with user-specified parameters (and equal weights for all points)
  sqpnp::SolverParameters params;
  params.omega_nullspace_method = sqpnp::OmegaNullspaceMethod::SVD;
  sqpnp::PnPSolver solver(points3d, points2d, std::vector<double>(n, 1.0), params);
```

## Cite as
If you use this code in your published work, please cite the following paper:<br><br>
@inproceedings{terzakis2020SQPnP,<br>
  title={A Consistently Fast and Globally Optimal Solution to the Perspective-n-Point Problem},<br>
  author={George Terzakis and Manolis Lourakis},<br>
  booktitle={European Conference on Computer Vision},<br>
  pages={478--494},<br>
  year={2020},<br>
  publisher={Springer International Publishing}<br>
}<br>
