//
// sqpnp.h
//
// George Terzakis (terzakig-at-hotmail-dot-com), September 2020
//
// Implementation of SQPnP as described in the paper:
//
// "A Consistently Fast and Globally Optimal Solution to the Perspective-n-Point Problem" by G. Terzakis and M. Lourakis
//  	 a) Paper: 	   http://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460.pdf 
//       b) Supplementary: https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460.pdf
  

#ifndef SQPnP_H__
#define SQPnP_H__

#include "types.h"
#include <vector>
#include <assert.h>
#include <iostream>

#include <Eigen/SparseCholesky> 
#include <Eigen/Sparse>

namespace sqpnp
{
  
  class PnPSolver
  {
    
  public:
    
    static const double SQRT3;
    
    bool IsValid() const { return flag_valid_; }
    const Eigen::Matrix<double, 9, 9>& Omega() const { return Omega_; }
    const Eigen::Matrix<double, 9, 9>& EigenVectors() const { return U_; }
    const Eigen::Vector<double, 9>& EigenValues() const { return s_; }
    const int NullSpaceDimension() const { return num_null_vectors_; }
    const int NumberOfSolutions() const { return num_solutions_; }
    const SQPSolution* const SolutionPtr(const int index) const 
    { 
      return index < 0 || index > num_solutions_-1 ? nullptr : &solutions_[index]; 
    }
    
    //
    // Constructor (initializes Omega and P and U, s, i.e. the decomposition of Omega)
    template <class Point3D, class Projection2D>
    inline PnPSolver(const std::vector<Point3D>& _3dpoints, 
		    const std::vector<Projection2D>& _projections, 
		    const SolverParameters& _parameters = SolverParameters()
		    ) : parameters_(_parameters)
    {
      if (_3dpoints.size() !=_projections.size() || _3dpoints.size() < 3 || _projections.size() < 3 )
      {
	flag_valid_= false;
	return;
      }
      
      flag_valid_ = true;
      num_null_vectors_ = -1; // set to -1 in case we never make it to the decomposition of Omega
      Omega_ = Eigen::Matrix<double, 9, 9>::Zero();
      const size_t n = _3dpoints.size();
      double sum_x = 0, 
	     sum_y = 0, 
	     sum_x2_plus_y2 = 0;
      double sum_X = 0,
	     sum_Y = 0,
	     sum_Z = 0;
	     
      Eigen::Matrix<double, 3, 9> QA = Eigen::Matrix<double, 3, 9>::Zero();  // Sum( Qi*Ai )
      
      for (size_t i = 0; i < n; i++)
      {
	points_.push_back( _3dpoints.at(i) );
	projections_.push_back( _projections.at(i) );
	
	double x = projections_.rbegin()->vector[0],
	       y = projections_.rbegin()->vector[1], 
	       sq_norm_m = projections_.rbegin()->vector.squaredNorm();
	sum_x += x;
	sum_y += y;
	sum_x2_plus_y2 += sq_norm_m;
	
	double X = points_.rbegin()->vector[0],
	       Y = points_.rbegin()->vector[1],
	       Z = points_.rbegin()->vector[2];
	sum_X += X;
	sum_Y += Y;
	sum_Z += Z;
	
	// Accumulate Omega by kronecker( Qi, Mi*Mi' ) = Q + A'*Qi*Ai. NOTE: Skipping block (3:5, 3:5) because its same as (0:2, 0:2)
	double X2 = X*X, XY = X*Y, XZ = X*Z, Y2 = Y*Y, YZ = Y*Z, Z2 = Z*Z;
	// a. Block (0:2, 0:2) populated by Mi*Mi'. NOTE: Only upper triangle
	Omega_(0, 0) += X2;
	Omega_(0, 1) += XY;
	Omega_(0, 2) += XZ;
	Omega_(1, 1) += Y2;
	Omega_(1, 2) += YZ;
	Omega_(2, 2) += Z2;
	
	// b. Block (0:2, 6:8) populated by -x*Mi*Mi'. NOTE: Only upper triangle
	Omega_(0, 6) += -x*X2; Omega_(0, 7) += -x*XY; Omega_(0, 8) += -x*XZ;
			       Omega_(1, 7) += -x*Y2; Omega_(1, 8) += -x*YZ;  
						      Omega_(2, 8) += -x*Z2;
	// c. Block (3:5, 6:8) populated by -y*Mi*Mi'. NOTE: Only upper triangle
	Omega_(3, 6) += -y*X2; Omega_(3, 7) += -y*XY; Omega_(3, 8) += -y*XZ;
			       Omega_(4, 7) += -y*Y2; Omega_(4, 8) += -y*YZ;  
						      Omega_(5, 8) += -y*Z2;
	// d. Block (6:8, 6:8) populated by (x^2+y^2)*Mi*Mi'. NOTE: Only upper triangle
	Omega_(6, 6) += sq_norm_m*X2; Omega_(6, 7) += sq_norm_m*XY; Omega_(6, 8) += sq_norm_m*XZ;
				      Omega_(7, 7) += sq_norm_m*Y2; Omega_(7, 8) += sq_norm_m*YZ;  
								    Omega_(8, 8) += sq_norm_m*Z2;
									
	// Accumulating Qi*Ai in QA
	QA(0, 0) += X; QA(0, 1) += Y; QA(0, 2) += Z; 	QA(0, 6) += -x*X; QA(0, 7) += -x*Y; QA(0, 8) += -x*Z;
	QA(1, 3) += X; QA(1, 4) += Y; QA(1, 5) += Z; 	QA(1, 6) += -y*X; QA(1, 7) += -y*Y; QA(1, 8) += -y*Z;
	
	QA(2, 0) += -x*X; QA(2, 1) += -x*Y; QA(2, 2) += -x*Z; 	QA(2, 3) += -y*X; QA(2, 4) += -y*Y; QA(2, 5) += -y*Z;
	QA(2, 6) += sq_norm_m*X; QA(2, 7) += sq_norm_m*Y; QA(2, 8) += sq_norm_m*Z; 
	
      }
      
      // Fill-in lower triangles of off-diagonal blocks (0:2, 6:8), (3:5, 6:8) and (6:8, 6:8)
      Omega_(1, 6) = Omega_(0, 7); Omega_(2, 6) = Omega_(0, 8); Omega_(2, 7) = Omega_(1, 8);
      Omega_(4, 6) = Omega_(3, 7); Omega_(5, 6) = Omega_(3, 8); Omega_(5, 7) = Omega_(4, 8);
      Omega_(7, 6) = Omega_(6, 7); Omega_(8, 6) = Omega_(6, 8); Omega_(8, 7) = Omega_(7, 8);
      
      // Fill-in upper triangle of block (3:5, 3:5)
      Omega_(3, 3) = Omega_(0, 0); Omega_(3, 4) = Omega_(0, 1); Omega_(3, 5) = Omega_(0, 2);
				   Omega_(4, 4) = Omega_(1, 1); Omega_(4, 5) = Omega_(1, 2);
								Omega_(5, 5) = Omega_(2, 2);
      // Fill lower triangle of Omega
      for (int r = 0; r < 9; r++)
      {
	for (int c = 0; c < r; c++)
	{
	  Omega_(r, c) = Omega_(c, r);
	}
      }
       
      // Q = Sum( Qi ) = Sum( [ 1, 0, -x; 0, 1, -y; -x, -y, x^2 + y^2] )
      Eigen::Matrix<double, 3, 3> Q;
      Q(0, 0) = n; Q(0, 1) = 0; Q(0, 2) = -sum_x;
      Q(1, 0) = 0; Q(1, 1) = n; Q(1, 2) = -sum_y;
      Q(2, 0) = -sum_x; Q(2, 1) = -sum_y; Q(2, 2) = sum_x2_plus_y2; 
      
      // Qinv = inv( Q ) = inv( Sum( Qi) )
      double inv_detQ = 1.0 / ( n*( n*sum_x2_plus_y2 - sum_y*sum_y - sum_x*sum_x ) );
      if ( inv_detQ < 1e-6 ) 
      {
	flag_valid_ = false;
	return;
      }
      
      Eigen::Matrix<double, 3, 3> Qinv;
      Qinv(0, 0) =  inv_detQ * ( n*sum_x2_plus_y2 - sum_y*sum_y ); 
      Qinv(0, 1) = Qinv(1, 0) = inv_detQ * sum_x*sum_y;
      Qinv(0, 2) = Qinv(2, 0) = inv_detQ * n*sum_x;
      Qinv(1, 1) = inv_detQ * ( n*sum_x2_plus_y2 - sum_x*sum_x );
      Qinv(1, 2) = Qinv(2, 1) = inv_detQ * n*sum_y;
      Qinv(2, 2) = inv_detQ * n*n;
      
      // Compute P = -inv( Sum(Qi) ) * Sum( Qi*Ai ) = -Qinv * QA
      P_ = -Qinv * QA;
      // Complete Omega (i.e., Omega = Sum(A'*Qi*A') + Sum(Qi*Ai)'*P = Sum(A'*Qi*A') + Sum(Qi*Ai)'*inv(Sum(Qi))*Sum( Qi*Ai) 
      Omega_ +=  QA.transpose()*P_;
      
      // Finally, decompose Omega
      Eigen::JacobiSVD<Eigen::Matrix<double, 9, 9>> svd(Omega_, Eigen::ComputeFullU);
        
      U_ = svd.matrixU();
      s_ = svd.singularValues();
      if ( s_[0] < 1e-7 ) 
      {
	flag_valid_ = false;
	return;
      }
      // Find dimension of null space
      while (s_[7 - num_null_vectors_] / s_[0] < _parameters.rank_tolerance) num_null_vectors_++;
      // Dimension of null space of Omega must be <= 6
      if (++num_null_vectors_ > 6) 
      {
	flag_valid_ = false;
      }
      
      // Point mean 
      double inv_n = 1.0 / n;
      point_mean_ << sum_X*inv_n, sum_Y*inv_n, sum_Z*inv_n;
    }
    
    // Solve the PnP
    bool Solve();
    
  private:
    std::vector<_Projection> projections_;
    std::vector<_Point> points_;
    
    SolverParameters parameters_;
    
    Eigen::Matrix<double, 9, 9> Omega_;
    
    Eigen::Vector<double, 9> s_;
    Eigen::Matrix<double, 9, 9> U_;
    Eigen::Matrix<double, 3, 9> P_;
    Eigen::Vector<double, 3> point_mean_; // For the positive depth test
    
    int num_null_vectors_;
    
    bool flag_valid_;
    
    SQPSolution solutions_[18];
    int num_solutions_;
    
    //
    // Run sequential quadratic programming with orthogonality constraints
    SQPSolution RunSQP(const Eigen::Vector<double, 9>& r0);
    
    //
    // Solve the SQP system efficiently.
    void SolveSQPSystem(const Eigen::Vector<double, 9>& r, Eigen::Vector<double, 9>& delta );
    
    // Handle a newly found solution and populate the list of solutions
    void HandleSolution(SQPSolution& solution, double& min_sq_error);
    
    //
    // Test cheirality for a given solution
    inline bool TestPositiveDepth(const SQPSolution& solution)
    {
      const auto& r = solution.r_hat;
      const auto& t = solution.t;
      const auto& M = point_mean_;
      return ( r[6]*M[0] + r[7]*M[1] + r[8]*M[2] + t[2] > 0 );
	
    }
    
    //
    // Determinant of 3x3 matrix stored in a vector in ROW-MAJOR fashion
    inline static double Determinant3x3(const Eigen::Vector<double, 9>& r)
    {
      return r[0]*r[4]*r[8] + r[1]*r[5]*r[6] + r[2]*r[3]*r[7] - r[6]*r[4]*r[2] - r[7]*r[5]*r[0] - r[8]*r[3]*r[1];
    }
    
    //
    // Invert a 3x3 symmetrix matrix (using low triangle values only)
    inline static bool InvertSymmetric3x3(const Eigen::Matrix<double, 3, 3> Q, 
					  Eigen::Matrix<double, 3, 3>& Qinv, 
					  const double& det_threshold = 1e-8
					 )
    {
       // 1. Get the elements of the matrix
       double a = Q(0, 0), b = Q(1, 0), c = Q(2, 0),
	      d = Q(1, 0), e = Q(1, 1), f = Q(2, 1),
	      g = Q(2, 0), h = Q(2, 1), i = Q(2, 2);
 
      // 2. Determinant
      double det = a*e*i + b*f*g + c*d*h - g*e*c - h*f*a - i*d*b;
      
      if ( fabs(det) < det_threshold ) return false;
      double invDet = 1.0 / det;
 
      // 3. Adjoint and inverse
      Qinv(0, 0) = ( e*i - f*h ) * invDet; Qinv(0, 1) = ( h*c - i*b ) * invDet;   Qinv(0, 2) = ( b*f - c*e ) * invDet;
      Qinv(1, 0) =  Qinv(0, 1);            Qinv(1, 1) = ( a*i - g*c ) * invDet;   Qinv(1, 2) = ( d*c - a*f ) * invDet;
      Qinv(2, 0) =  Qinv(0, 2);            Qinv(2, 1) = Qinv(1, 2);               Qinv(2, 2) = ( a*e - d*b ) * invDet;
    }
    
    // Simple SVD - based nearest rotation matrix. Argument should be a ROW-MAJOR matrix representation.
    // Returns a ROW-MAJOR vector representation of the nearest rotation matrix.
    // NOTE: This could be improved by adding the OLAE or FOAMOA methods ( see: http://users.ics.forth.gr/~lourakis/publ/2018_iros.pdf ) 
    inline static void NearestRotationMatrix(const Eigen::Vector<double, 9>& e, Eigen::Vector<double, 9>& r)
    {
      const Eigen::Matrix<double, 3, 3> E = e.reshaped(3, 3); 
      Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3>> svd( E, Eigen::ComputeFullU | Eigen::ComputeFullV );
      double detUV = Determinant3x3(svd.matrixU().reshaped(9, 1)) * Determinant3x3(svd.matrixV().reshaped(9, 1)); // Row/col major doesn't play a role here...
      // so we return back a row-major vector representation of the orthogonal matrix
      r = ( svd.matrixU() * Eigen::Vector<double, 3>({1, 1, detUV}).asDiagonal() * svd.matrixV().transpose() ).reshaped(9, 1);
    }
    
    //
    // Produce a distance from being orthogonal for a random 3x3 matrix
    // Matrix is provided as a vector
    inline static double OrthogonalityError(const Eigen::Vector<double, 9>& a)
    {
      double sq_norm_a1 = a[0]*a[0] + a[1]*a[1] + a[2]*a[2],
             sq_norm_a2 = a[3]*a[3] + a[4]*a[4] + a[5]*a[5],
	     sq_norm_a3 = a[6]*a[6] + a[7]*a[7] + a[8]*a[8];
      double dot_a1a2 = a[0]*a[3] + a[1]*a[4] + a[2]*a[5],
	     dot_a1a3 = a[0]*a[6] + a[1]*a[7] + a[2]*a[8],
	     dot_a2a3 = a[3]*a[6] + a[4]*a[7] + a[5]*a[8];
	     
      return (sq_norm_a1 - 1)*(sq_norm_a1 - 1) + (sq_norm_a2 - 1)*(sq_norm_a2 - 1) + (sq_norm_a3 - 1)*(sq_norm_a3 - 1) + 
	      2*( dot_a1a2*dot_a1a2 + dot_a1a3*dot_a1a3 + dot_a2a3*dot_a2a3 );
    }
    
    //
    // Compute the 3D null space and 6D normal space of the constraint Jacobian at a 9D vector r (not necessarilly a rotation-yet it should be rank-3)
    static void RowAndNullSpace(const Eigen::Vector<double, 9>& r, 
				      Eigen::Matrix<double, 9, 6>& Q, 
				      Eigen::Matrix<double, 9, 3>& N,
				      Eigen::Matrix<double, 6, 6>& K,
				const double& norm_threhsold = 0.1 );
    
    
  };
  

}

#endif
