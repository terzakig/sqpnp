 //
// Types.h
//
// George Terzakis, September 2020
//
// Implementation of SQPnP as described in the paper:
//
// "A Consistently Fast and Globally Optimal Solution to the Perspective-n-Point Problem" by G. Terzakis and M. Lourakis
//  	 a) Paper: 	   http://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460.pdf 
//       b) Supplementary: https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460.pdf

#ifndef _SQPNPTYPES__H_
#define _SQPNPTYPES__H_
 
  
#ifdef HAVE_OPENCV

#include <opencv2/core.hpp>

#endif

#include <iostream>
#include <Eigen/Dense>
 
namespace sqpnp
{
 
  // Euclidean projection
  struct _Projection
  {
    Eigen::Matrix<double, 2, 1> vector;
    
    inline _Projection() : vector(Eigen::Matrix<double, 2, 1>(0, 0)) {}
    
    template<typename P>
    inline _Projection(const P& _x, const P& _y) : vector(Eigen::Matrix<double, 2, 1>(_x, _y)) {}
    
    template <class Proj>
    inline _Projection(const Proj& _projection) { *this = _projection; }
    
    #ifdef HAVE_OPENCV
    template <typename P>
    _Projection& operator =(const cv::Point_<P>& _projection) 
    {
	vector[0] = _projection.x;
	vector[1] = _projection.y;
	return *this;
    }
    template <typename P>
    _Projection& operator =(const cv::Vec<P, 2>& _projection) 
    {
	vector[0] = _projection[0];
	vector[1] = _projection[1];
	return *this;
    }
    #endif
    
    template <typename P>
    _Projection& operator =(const Eigen::Matrix<P, 2, 1>& _projection) 
    {
	vector = _projection;
	return *this;
    }

  };
 
  // 3D point
  struct _Point
  {
    Eigen::Matrix<double, 3, 1> vector;
    
    inline _Point() : vector(Eigen::Matrix<double, 3, 1>(0, 0, 0)) {}
    
    template<typename P>
    inline _Point( const P& _x, const P& _y, const P& _z ) : vector(Eigen::Matrix<double, 3, 1>(_x, _y, _z)) {}
    
    template <class Pt>
    inline _Point(const Pt& _point) { *this = _point; }
    
    #ifdef HAVE_OPENCV
    template <typename P>
    _Point& operator =(const cv::Point3_<P>& _point) 
    {
	vector[0] = _point.x;
	vector[1] = _point.y;
	vector[2] = _point.z;
	return *this;
    }
    template <typename P>
    _Point& operator =(const cv::Vec<P, 3>& _point) 
    {
	vector[0] = _point[0];
	vector[1] = _point[1];
	vector[2] = _point[2];
	return *this;
    }
    #endif
    
    template <typename P>
    _Point& operator =(const Eigen::Matrix<P, 3, 1>& _point) 
    {
	vector = _point;
	return *this;
    }

  };
  
  struct SolverParameters
  {
    static const double DEFAULT_RANK_TOLERANCE;
    static const double DEFAULT_SQP_SQUARED_TOLERANCE;
    static const double DEFAULT_SQP_DET_THRESHOLD;
    static const int DEFAULT_SQP_MAX_ITERATION = 15;
    static const double DEFAULT_ORTHOGONALITY_SQUARED_ERROR_THRESHOLD;
    static const double DEFAULT_EQUAL_VECTORS_SQUARED_DIFF;
    static const double DEFAULT_EQUAL_SQUARED_ERRORS_DIFF;
    static const double DEFAULT_POINT_VARIANCE_THRESHOLD;
    
    double rank_tolerance;
    double sqp_squared_tolerance;
    double sqp_det_threshold;
    int sqp_max_iteration;
    double orthogonality_squared_error_threshold;  
    double equal_vectors_squared_diff;
    double equal_squared_errors_diff;
    double point_variance_threshold;
    
    inline SolverParameters(const double& _rank_tolerance = DEFAULT_RANK_TOLERANCE, 
			    const double& _sqp_squared_tolerance = DEFAULT_SQP_SQUARED_TOLERANCE,
			    const double& _sqp_det_threshold = DEFAULT_SQP_DET_THRESHOLD,
			    const int _sqp_max_iteration = DEFAULT_SQP_MAX_ITERATION, 
			    const double& _orthogonality_squared_error_threshold = DEFAULT_ORTHOGONALITY_SQUARED_ERROR_THRESHOLD,
		            const double& _equal_vectors_squared_diff = DEFAULT_EQUAL_VECTORS_SQUARED_DIFF,
			    const double& _equal_squared_errors_diff = DEFAULT_EQUAL_SQUARED_ERRORS_DIFF,
			    const double& _point_variance_threshold = DEFAULT_POINT_VARIANCE_THRESHOLD
			    ) : rank_tolerance(_rank_tolerance),
				sqp_squared_tolerance(_sqp_squared_tolerance),
				sqp_det_threshold(_sqp_det_threshold),
				sqp_max_iteration(_sqp_max_iteration),
				orthogonality_squared_error_threshold(_orthogonality_squared_error_threshold),
				equal_vectors_squared_diff(_equal_vectors_squared_diff),
				equal_squared_errors_diff(_equal_squared_errors_diff),
				point_variance_threshold(_point_variance_threshold)
  {}
  
  };
  
  struct SQPSolution
  {
    Eigen::Matrix<double, 9, 1> r; // Actual matrix upon convergence
    Eigen::Matrix<double, 9, 1> r_hat; // "Clean" (nearest) rotation matrix 
    Eigen::Matrix<double, 3, 1> t;
    int num_iterations;
    double sq_error;
    
  };
  
  inline std::ostream& operator << (std::ostream& os, const SQPSolution& solution)
  {
    return os << "\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n" <<
	         "R: [ " << solution.r_hat[0] << ", " << solution.r_hat[1] << ", " << solution.r_hat[2] << ";\n" << 
		 "     " << solution.r_hat[3] << ", " << solution.r_hat[4] << ", " << solution.r_hat[5] << ";\n" <<
		 "     " << solution.r_hat[6] << ", " << solution.r_hat[7] << ", " << solution.r_hat[8] << " ]\n" <<
	         "t: " << "[ " << solution.t[0] << "; " << solution.t[1] << "; " << solution.t[2] << " ]\n" << 
	         "r: " << "[ " << solution.r_hat[0] << "; " << solution.r_hat[1] << "; " << solution.r_hat[2] << "; " 
			  << solution.r_hat[3] << "; " << solution.r_hat[4] << "; " << solution.r_hat[5] << "; " 
  			  << solution.r_hat[6] << "; " << solution.r_hat[7] << "; " << solution.r_hat[8] << " ]\n" <<
		 "squared error: " << solution.sq_error << "\n" <<
		 "number of SQP iterations : " << solution.num_iterations << "\n" <<
		 "-------------------------------------------------------------------------------------------------\n";
  }
	
  

}

#endif
