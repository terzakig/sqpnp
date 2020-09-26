# SQPnP 
C++ Implementation of the SQPnP algorithm. 

The algorithm is the generic PnP solver described in the paper ["A Consistently Fast and Globally Optimal Solution to the Perspective-n-Point Problem" by G. Terzakis and M. Lourakis](http://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460.pdf). For more intuition, refer to the supplementary material [here](https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460-supp.pdf).

## Required libraries
SQPnP requires Eigen to build. Besides the SVD, the use of the library is confined to matrix addition, transposition and multiplication. Choosing Eigen was motivated by its increasing popularity and lightweight character. There are two examles in this repository, a) one that uses OpenCV, just for the sake of demonstrating the initialization of SQPnP with ``cv::Point_<>`` and ``cv::Point3_<>`` structures, and b) a second example, in which data points are copied from plain 2D arrays. Build will proceed with either one, depending on whether OpenCV is found or not.

Build
-----

Create a ``build`` directory in the root of the cloned repository and run ``cmake``:

``mkdir build``

``cd build``

``cmake ..``

or, for a *release* build,

``cmake .. -DCMAKE_BUILD_TYPE=Release``

The latter will allow for more accurate timing of average execution time. Finally build everything:

``make``

To run the example, once in the ``build`` directory,

``./example/sqpnp_example``  
