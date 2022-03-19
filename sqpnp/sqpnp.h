//
// sqpnp.h
//
// George Terzakis (terzakig-at-hotmail-dot-com), September 2020
// Nearest orthogonal approximation code (C) 2019 Manolis Lourakis
//
// Implementation of SQPnP as described in the paper:
//
// "A Consistently Fast and Globally Optimal Solution to the Perspective-n-Point Problem" by G. Terzakis and M. Lourakis
//     a) Paper:         https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460.pdf 
//     b) Supplementary: https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460-supp.pdf
  

#ifndef SQPnP_H__
#define SQPnP_H__

#include "types.h"
#include <vector>
#include <assert.h>
#include <iostream>

namespace sqpnp
{
  
  class PnPSolver
  {
    
  public:
    
    static const double SQRT3;
    
    bool IsValid() const { return flag_valid_; }
    const Eigen::Matrix<double, 9, 9>& Omega() const { return Omega_; }
    const Eigen::Matrix<double, 9, 9>& EigenVectors() const { return U_; }
    const Eigen::Matrix<double, 9, 1>& EigenValues() const { return s_; }
    const int NullSpaceDimension() const { return num_null_vectors_; }
    const int NumberOfSolutions() const { return num_solutions_; }
    const SQPSolution* const SolutionPtr(const int index) const 
    { 
      return index < 0 || index >= num_solutions_ ? nullptr : &solutions_[index]; 
    }
    
    //
    // Return average reprojection errors
    inline std::vector<double> AverageSquaredProjectionErrors() const
    {
      std::vector<double> avg_errors;
      avg_errors.reserve(num_solutions_);
      for (int i = 0; i < num_solutions_; i++)
      {
          avg_errors.emplace_back( AverageSquaredProjectionError(i) );
      }
      return avg_errors;
    }
    const std::vector<double>& Weights() const { return weights_; }

    //
    // Constructor (initializes Omega, P and U, s, i.e. the decomposition of Omega)
    template <class Point3D, class Projection2D, typename Pw = double>
    inline PnPSolver(const std::vector<Point3D>& _3dpoints, 
		    const std::vector<Projection2D>& _projections, 
		    const std::vector<Pw>& _weights = std::vector<Pw>(),
		    const SolverParameters& _parameters = SolverParameters()) : parameters_(_parameters)
    {
      const size_t n = _3dpoints.size();

      if (n !=_projections.size() || n < 3 )
      {
        flag_valid_ = false;
        return;
      }
      
      if (!_weights.empty())
      {
        if (n != _weights.size() ) {
          flag_valid_ = false;
          return;
        }

        weights_ = _weights;
      }
      else
      {
        weights_.resize(n, 1.0);
      }
      
      flag_valid_ = true;
      points_.reserve(n);
      projections_.reserve(n);
      num_null_vectors_ = -1; // set to -1 in case we never make it to the decomposition of Omega
      Omega_ = Eigen::Matrix<double, 9, 9>::Zero();
      double sum_wx = 0.0,
        sum_wy = 0.0,
        sum_wx2_plus_wy2 = 0.0,
        sum_w = 0.0;
      
      double sum_wX = 0.0,
        sum_wY = 0.0,
        sum_wZ = 0.0;
	     
      Eigen::Matrix<double, 3, 9> QA = Eigen::Matrix<double, 3, 9>::Zero();  // Sum( Qi*Ai )
      
      for (size_t i = 0; i < n; i++)
      {
        const double w = weights_[i];
        points_.emplace_back( _3dpoints[i] );
        projections_.emplace_back( _projections[i] );

        if (w == 0.0) continue;
	
        const auto& proj_ = projections_[i].vector;
        const double wx = proj_[0] * w,
               wy = proj_[1] * w, 
        wsq_norm_m = w*proj_.squaredNorm();
        sum_wx += wx;
        sum_wy += wy;
        sum_wx2_plus_wy2 += wsq_norm_m;
        sum_w += w;

        const auto& pt_ = points_[i].vector;
        double X = pt_[0],
            Y = pt_[1],
            Z = pt_[2];
        const double wX = w * X,  wY = w * Y,  wZ = w * Z;
        sum_wX += wX;
        sum_wY += wY;
        sum_wZ += wZ;
	
        // Accumulate Omega by kronecker( Qi, Mi*Mi' ) = A'*Qi*Ai. NOTE: Skipping block (3:5, 3:5) because it's same as (0:2, 0:2)
        const double X2 = X*X, XY = X*Y, XZ = X*Z, Y2 = Y*Y, YZ = Y*Z, Z2 = Z*Z;
        
        // a. Block (0:2, 0:2) populated by Mi*Mi'. NOTE: Only upper triangle
        Omega_(0, 0) += w*X2;
        Omega_(0, 1) += w*XY;
        Omega_(0, 2) += w*XZ;
        Omega_(1, 1) += w*Y2;
        Omega_(1, 2) += w*YZ;
        Omega_(2, 2) += w*Z2;
        
        // b. Block (0:2, 6:8) populated by -x*Mi*Mi'. NOTE: Only upper triangle
        Omega_(0, 6) += -wx*X2; Omega_(0, 7) += -wx*XY; Omega_(0, 8) += -wx*XZ;
                                Omega_(1, 7) += -wx*Y2; Omega_(1, 8) += -wx*YZ;
                                                        Omega_(2, 8) += -wx*Z2;

        // c. Block (3:5, 6:8) populated by -y*Mi*Mi'. NOTE: Only upper triangle
        Omega_(3, 6) += -wy*X2; Omega_(3, 7) += -wy*XY; Omega_(3, 8) += -wy*XZ;
                                Omega_(4, 7) += -wy*Y2; Omega_(4, 8) += -wy*YZ;
                                                        Omega_(5, 8) += -wy*Z2;
                                                        
        // d. Block (6:8, 6:8) populated by (x^2+y^2)*Mi*Mi'. NOTE: Only upper triangle
        Omega_(6, 6) += wsq_norm_m*X2; Omega_(6, 7) += wsq_norm_m*XY; Omega_(6, 8) += wsq_norm_m*XZ;
                                       Omega_(7, 7) += wsq_norm_m*Y2; Omega_(7, 8) += wsq_norm_m*YZ;
                                                                      Omega_(8, 8) += wsq_norm_m*Z2;

        // Accumulating Qi*Ai in QA.
        // Note that certain pairs of elements are equal, so we save some operations by filling them outside the loop
        QA(0, 0) += wX; QA(0, 1) += wY; QA(0, 2) += wZ; 	QA(0, 6) += -wx*X; QA(0, 7) += -wx*Y; QA(0, 8) += -wx*Z;
        //QA(1, 3) += wX; QA(1, 4) += wY; QA(1, 5) += wZ;
                                                          QA(1, 6) += -wy*X; QA(1, 7) += -wy*Y; QA(1, 8) += -wy*Z;
        
        //QA(2, 0) += -wx*X; QA(2, 1) += -wx*Y; QA(2, 2) += -wx*Z; 	QA(2, 3) += -wy*X; QA(2, 4) += -wy*Y; QA(2, 5) += -wy*Z;
        QA(2, 6) += wsq_norm_m*X; QA(2, 7) += wsq_norm_m*Y; QA(2, 8) += wsq_norm_m*Z; 
      }
      // Complete QA
      QA(1, 3) = QA(0, 0); QA(1, 4) = QA(0, 1); QA(1, 5) = QA(0, 2);
      QA(2, 0) = QA(0, 6); QA(2, 1) = QA(0, 7); QA(2, 2) = QA(0, 8);
      QA(2, 3) = QA(1, 6); QA(2, 4) = QA(1, 7); QA(2, 5) = QA(1, 8);
      
      // Fill-in lower triangles of off-diagonal blocks (0:2, 6:8), (3:5, 6:8) and (6:8, 6:8)
      Omega_(1, 6) = Omega_(0, 7); Omega_(2, 6) = Omega_(0, 8); Omega_(2, 7) = Omega_(1, 8);
      Omega_(4, 6) = Omega_(3, 7); Omega_(5, 6) = Omega_(3, 8); Omega_(5, 7) = Omega_(4, 8);
      Omega_(7, 6) = Omega_(6, 7); Omega_(8, 6) = Omega_(6, 8); Omega_(8, 7) = Omega_(7, 8);
      
      // Fill-in upper triangle of block (3:5, 3:5)
      Omega_(3, 3) = Omega_(0, 0); Omega_(3, 4) = Omega_(0, 1); Omega_(3, 5) = Omega_(0, 2);
                                   Omega_(4, 4) = Omega_(1, 1); Omega_(4, 5) = Omega_(1, 2);
                                                                Omega_(5, 5) = Omega_(2, 2);

      // Fill lower triangle of Omega; elements (7, 6), (8, 6) & (8, 7) have already been assigned above
      Omega_(1, 0) = Omega_(0, 1);
      Omega_(2, 0) = Omega_(0, 2); Omega_(2, 1) = Omega_(1, 2);
      Omega_(3, 0) = Omega_(0, 3); Omega_(3, 1) = Omega_(1, 3); Omega_(3, 2) = Omega_(2, 3);
      Omega_(4, 0) = Omega_(0, 4); Omega_(4, 1) = Omega_(1, 4); Omega_(4, 2) = Omega_(2, 4); Omega_(4, 3) = Omega_(3, 4);
      Omega_(5, 0) = Omega_(0, 5); Omega_(5, 1) = Omega_(1, 5); Omega_(5, 2) = Omega_(2, 5); Omega_(5, 3) = Omega_(3, 5); Omega_(5, 4) = Omega_(4, 5);
      Omega_(6, 0) = Omega_(0, 6); Omega_(6, 1) = Omega_(1, 6); Omega_(6, 2) = Omega_(2, 6); Omega_(6, 3) = Omega_(3, 6); Omega_(6, 4) = Omega_(4, 6); Omega_(6, 5) = Omega_(5, 6);
      Omega_(7, 0) = Omega_(0, 7); Omega_(7, 1) = Omega_(1, 7); Omega_(7, 2) = Omega_(2, 7); Omega_(7, 3) = Omega_(3, 7); Omega_(7, 4) = Omega_(4, 7); Omega_(7, 5) = Omega_(5, 7);
      Omega_(8, 0) = Omega_(0, 8); Omega_(8, 1) = Omega_(1, 8); Omega_(8, 2) = Omega_(2, 8); Omega_(8, 3) = Omega_(3, 8); Omega_(8, 4) = Omega_(4, 8); Omega_(8, 5) = Omega_(5, 8);

      // Q = Sum( wi*Qi ) = Sum( [ wi, 0, -wi*xi; 0, 1, -wi*yi; -wi*xi, -wi*yi, wi*(xi^2 + yi^2)] )
      Eigen::Matrix<double, 3, 3> Q;
      Q(0, 0) = sum_w;     Q(0, 1) = 0.0;       Q(0, 2) = -sum_wx;
      Q(1, 0) = 0.0;       Q(1, 1) = sum_w;     Q(1, 2) = -sum_wy;
      Q(2, 0) = -sum_wx;   Q(2, 1) = -sum_wy;   Q(2, 2) = sum_wx2_plus_wy2; 
      
      // Qinv = inv( Q ) = inv( Sum( Qi) )
      Eigen::Matrix<double, 3, 3> Qinv;
      InvertSymmetric3x3(Q, Qinv);
      
      // Compute P = -inv( Sum(wi*Qi) ) * Sum( wi*Qi*Ai ) = -Qinv * QA
      P_ = -Qinv * QA;
      // Complete Omega (i.e., Omega = Sum(A'*Qi*A') + Sum(Qi*Ai)'*P = Sum(A'*Qi*A') + Sum(Qi*Ai)'*inv(Sum(Qi))*Sum( Qi*Ai) 
      Omega_ +=  QA.transpose()*P_;
      
      // Finally, decompose Omega with the chosen method
      if ( parameters_.omega_nullspace_method == OmegaNullspaceMethod::RRQR )
      {
        // Rank revealing QR nullspace computation. This is slightly less accurate compared to SVD but x2 faster
        Eigen::FullPivHouseholderQR<Eigen::Matrix<double, 9, 9> > rrqr(Omega_);
        U_ = rrqr.matrixQ();

        Eigen::Matrix<double, 9, 9> R = rrqr.matrixQR().template triangularView<Eigen::Upper>();
        s_ = R.diagonal().array().abs();
      }
      else // if ( parameters_.omega_nullspace_method == OmegaNullspaceMethod::SVD )
      {
        // SVD-based nullspace computation. This is the most accurate but slowest option
        Eigen::JacobiSVD<Eigen::Matrix<double, 9, 9>> svd(Omega_, Eigen::ComputeFullU);
        U_ = svd.matrixU();
        s_ = svd.singularValues();
      }

      // Find dimension of null space; the check guards against overly large rank_tolerance
      while (7 - num_null_vectors_ >= 0 && s_[7 - num_null_vectors_] < _parameters.rank_tolerance) num_null_vectors_++;
      //while (s_[7 - num_null_vectors_] < _parameters.rank_tolerance) num_null_vectors_++;

      // Dimension of null space of Omega must be <= 6
      if (++num_null_vectors_ > 6) 
      {
          flag_valid_ = false;
      }
      
      // 3D point weighted mean (for quick cheirality checks)
      const double inv_sum_w = 1.0 / sum_w;
      point_mean_ << sum_wX*inv_sum_w, sum_wY*inv_sum_w, sum_wZ*inv_sum_w;
       
      // Assign nearest rotation method
      if ( parameters_.nearest_rotation_method == NearestRotationMethod::FOAM )
      {
          NearestRotationMatrix = NearestRotationMatrix_FOAM;
      }
      else // if ( parameters_.nearest_rotation_method == NearestRotationMethod::SVD )
      {
          NearestRotationMatrix = NearestRotationMatrix_SVD;
      }
    }
    
    // Solve the PnP
    bool Solve();
    void H(int arg1, int arg2);
    
  private:
    std::vector<_Projection> projections_;
    std::vector<_Point> points_;
    std::vector<double> weights_;
    SolverParameters parameters_;
    
    Eigen::Matrix<double, 9, 9> Omega_;
    
    Eigen::Matrix<double, 9, 1> s_;
    Eigen::Matrix<double, 9, 9> U_;
    Eigen::Matrix<double, 3, 9> P_;
    Eigen::Matrix<double, 3, 1> point_mean_; // For the positive depth test
    
    int num_null_vectors_;
    
    bool flag_valid_;
    
    SQPSolution solutions_[18];
    int num_solutions_;
    
    // Nearest rotation matrix function. By default, the FOAM method
    std::function<void(const Eigen::Matrix<double, 9, 1>&, Eigen::Matrix<double, 9, 1>&)> NearestRotationMatrix;
    
    //
    // Run sequential quadratic programming with orthogonality constraints
    SQPSolution RunSQP(const Eigen::Matrix<double, 9, 1>& r0);
    
    //
    // Solve the SQP system efficiently.
    void SolveSQPSystem(const Eigen::Matrix<double, 9, 1>& r, Eigen::Matrix<double, 9, 1>& delta );
    
    // Handle a newly found solution and populate the list of solutions
    void HandleSolution(SQPSolution& solution, double& min_sq_error);
	
    //
    // Average squared projection error of a given solution
    inline double AverageSquaredProjectionError(const int index) const
    {
      double avg = 0.0;
      const auto& r = solutions_[index].r_hat;
      const auto& t = solutions_[index].t;
      
      for (size_t i = 0; i < points_.size(); i++)
      {
	const auto& M = points_[i].vector;
	const double Xc     =         r[0]*M[0] + r[1]*M[1] + r[2]*M[2] + t[0], 
	       Yc     =         r[3]*M[0] + r[4]*M[1] + r[5]*M[2] + t[1], 
	       inv_Zc = 1.0 / ( r[6]*M[0] + r[7]*M[1] + r[8]*M[2] + t[2] );

	const auto& m = projections_[i].vector;
	const double dx = Xc*inv_Zc - m[0];
	const double dy = Yc*inv_Zc - m[1];
	avg += dx*dx + dy*dy;
      }
      
      return avg / points_.size();
    }    

    //
    // Test cheirality on the mean point for a given solution
    inline bool TestPositiveDepth(const SQPSolution& solution)
    {
      const auto& r = solution.r_hat;
      const auto& t = solution.t;
      const auto& M = point_mean_;
      return ( r[6]*M[0] + r[7]*M[1] + r[8]*M[2] + t[2] > 0 );
    }

    //
    // Test cheirality on the majority of points for a given solution
    inline bool TestPositiveMajorityDepths(const SQPSolution& solution)
    {
      const auto& r = solution.r_hat;
      const auto& t = solution.t;
      int npos = 0, nneg = 0;

      for (size_t i = 0; i < points_.size(); i++) {
        if ( weights_[i] == 0.0) continue;
        const auto& M = points_[i].vector;
        if ( r[6]*M[0] + r[7]*M[1] + r[8]*M[2] + t[2] > 0 ) ++npos;
        else ++nneg;
      }

      return npos >= nneg;
    }
    
    //
    // Determinant of 3x3 matrix stored as a 9x1 vector in *row-major* order
    inline static double Determinant9x1(const Eigen::Matrix<double, 9, 1>& r)
    {
      return r[0]*r[4]*r[8] + r[1]*r[5]*r[6] + r[2]*r[3]*r[7] - r[6]*r[4]*r[2] - r[7]*r[5]*r[0] - r[8]*r[3]*r[1];
    }

    // Determinant of 3x3 matrix
    inline static double Determinant3x3(const Eigen::Matrix<double, 3, 3>& M)
    {
      return M(0, 0)*( M(1, 1)*M(2, 2)  - M(1, 2)*M(2, 1) ) - M(0, 1)*( M(1, 0)*M(2, 2)  - M(1, 2)*M(2, 0) ) + M(0, 2)*( M(1, 0)*M(2, 1)  - M(1, 1)*M(2, 0) );
    }
    
    //
    // Invert a 3x3 symmetric matrix (using low triangle values only)
    inline static bool InvertSymmetric3x3(const Eigen::Matrix<double, 3, 3> Q, 
					  Eigen::Matrix<double, 3, 3>& Qinv, 
					  const double& det_threshold = 1e-10
					 )
    {
      // 1. Get the elements of the matrix
      double a = Q(0, 0),
             b = Q(1, 0), d = Q(1, 1),
             c = Q(2, 0), e = Q(2, 1), f = Q(2, 2);
 
      // 2. Determinant
      double t2, t4, t7, t9, t12;
      t2 = e*e;
      t4 = a*d;
      t7 = b*b;
      t9 = b*c;
      t12 = c*c;
      double det = -t4*f+a*t2+t7*f-2.0*t9*e+t12*d;
      
      if ( fabs(det) < det_threshold ) { Qinv=Q.completeOrthogonalDecomposition().pseudoInverse(); return false; } // fall back to pseudoinverse
 
      // 3. Inverse
      double t15, t20, t24, t30;
      t15 = 1.0/det;
      t20 = (-b*f+c*e)*t15;
      t24 = (b*e-c*d)*t15;
      t30 = (a*e-t9)*t15;
      Qinv(0, 0) = (-d*f+t2)*t15;
      Qinv(0, 1) = Qinv(1, 0) = -t20;
      Qinv(0, 2) = Qinv(2, 0) = -t24;
      Qinv(1, 1) = -(a*f-t12)*t15;
      Qinv(1, 2) = Qinv(2, 1) = t30;
      Qinv(2, 2) = -(t4-t7)*t15;

      return true;
    }
    

    // Simple SVD - based nearest rotation matrix. Argument should be a *row-major* matrix representation.
    // Returns a row-major vector representation of the nearest rotation matrix.
    inline static void NearestRotationMatrix_SVD(const Eigen::Matrix<double, 9, 1>& e, Eigen::Matrix<double, 9, 1>& r)
    {
      //const Eigen::Matrix<double, 3, 3> E = e.reshaped(3, 3); 
      const Eigen::Map<const Eigen::Matrix<double, 3, 3>> E(e.data(), 3, 3);
      Eigen::JacobiSVD<Eigen::Matrix<double, 3, 3>> svd( E, Eigen::ComputeFullU | Eigen::ComputeFullV );
      double detUV = Determinant3x3(svd.matrixU()) * Determinant3x3(svd.matrixV());
      // so we return back a row-major vector representation of the orthogonal matrix
      Eigen::Matrix<double, 3, 3> R = svd.matrixU() * Eigen::Matrix<double, 3, 1>({1, 1, detUV}).asDiagonal() * svd.matrixV().transpose();
      r = Eigen::Map<Eigen::Matrix<double, 9, 1>>(R.data(), 9, 1);
    }

    // Faster nearest rotation computation based on FOAM. See M. Lourakis: "An Efficient Solution to Absolute Orientation", ICPR 2016
    // and M. Lourakis, G. Terzakis: "Efficient Absolute Orientation Revisited", IROS 2018.
    /* Solve the nearest orthogonal approximation problem
     * i.e., given B, find R minimizing ||R-B||_F
     *
     * The computation borrows from Markley's FOAM algorithm
     * "Attitude Determination Using Vector Observations: A Fast Optimal Matrix Algorithm", J. Astronaut. Sci. 1993.
     *
     *  Copyright (C) 2019 Manolis Lourakis (lourakis **at** ics forth gr)
     *  Institute of Computer Science, Foundation for Research & Technology - Hellas
     *  Heraklion, Crete, Greece.
     */
    inline static void NearestRotationMatrix_FOAM(const Eigen::Matrix<double, 9, 1>& e, Eigen::Matrix<double, 9, 1>& r)
    {
    int i;
    const double *B=e.data();
    double l, lprev, detB, Bsq, adjBsq, adjB[9];

      // det(B)
      detB=B[0]*B[4]*B[8] - B[0]*B[5]*B[7] - B[1]*B[3]*B[8] + B[2]*B[3]*B[7] + B[1]*B[6]*B[5] - B[2]*B[6]*B[4];
      if(fabs(detB)<1E-04){ // singular, let SVD handle it
        NearestRotationMatrix_SVD(e, r);
        return;
      }

      // B's adjoint
      adjB[0]=B[4]*B[8] - B[5]*B[7]; adjB[1]=B[2]*B[7] - B[1]*B[8]; adjB[2]=B[1]*B[5] - B[2]*B[4];
      adjB[3]=B[5]*B[6] - B[3]*B[8]; adjB[4]=B[0]*B[8] - B[2]*B[6]; adjB[5]=B[2]*B[3] - B[0]*B[5];
      adjB[6]=B[3]*B[7] - B[4]*B[6]; adjB[7]=B[1]*B[6] - B[0]*B[7]; adjB[8]=B[0]*B[4] - B[1]*B[3];

      // ||B||^2, ||adj(B)||^2
      Bsq=B[0]*B[0]+B[1]*B[1]+B[2]*B[2] + B[3]*B[3]+B[4]*B[4]+B[5]*B[5] + B[6]*B[6]+B[7]*B[7]+B[8]*B[8];
      adjBsq=adjB[0]*adjB[0]+adjB[1]*adjB[1]+adjB[2]*adjB[2] + adjB[3]*adjB[3]+adjB[4]*adjB[4]+adjB[5]*adjB[5] + adjB[6]*adjB[6]+adjB[7]*adjB[7]+adjB[8]*adjB[8];

      // compute l_max with Newton-Raphson from FOAM's characteristic polynomial, i.e. eq.(23) - (26)
      l=0.5*(Bsq + 3.0); // 1/2*(trace(B*B') + trace(eye(3)))
      if(detB<0.0) l=-l; // trB & detB have opposite signs!
      for(i=15, lprev=0.0; fabs(l-lprev)>1E-12*fabs(lprev) && i>0; --i){
        double tmp, p, pp;

        tmp=(l*l-Bsq);
        p=(tmp*tmp - 8.0*l*detB - 4.0*adjBsq);
        pp=8.0*(0.5*tmp*l - detB);

        lprev=l;
        l-=p/pp;
      }

      // the rotation matrix equals ((l^2 + Bsq)*B + 2*l*adj(B') - 2*B*B'*B) / (l*(l*l-Bsq) - 2*det(B)), i.e. eq.(14) using (18), (19)
      {
      // compute (l^2 + Bsq)*B
      double tmp[9], BBt[9], denom;
      const double a=l*l + Bsq;

      // BBt=B*B'
      BBt[0]=B[0]*B[0] + B[1]*B[1] + B[2]*B[2];
      BBt[1]=B[0]*B[3] + B[1]*B[4] + B[2]*B[5];
      BBt[2]=B[0]*B[6] + B[1]*B[7] + B[2]*B[8];

      BBt[3]=BBt[1];
      BBt[4]=B[3]*B[3] + B[4]*B[4] + B[5]*B[5];
      BBt[5]=B[3]*B[6] + B[4]*B[7] + B[5]*B[8];

      BBt[6]=BBt[2];
      BBt[7]=BBt[5];
      BBt[8]=B[6]*B[6] + B[7]*B[7] + B[8]*B[8];

      // tmp=BBt*B
      tmp[0]=BBt[0]*B[0] + BBt[1]*B[3] + BBt[2]*B[6];
      tmp[1]=BBt[0]*B[1] + BBt[1]*B[4] + BBt[2]*B[7];
      tmp[2]=BBt[0]*B[2] + BBt[1]*B[5] + BBt[2]*B[8];

      tmp[3]=BBt[3]*B[0] + BBt[4]*B[3] + BBt[5]*B[6];
      tmp[4]=BBt[3]*B[1] + BBt[4]*B[4] + BBt[5]*B[7];
      tmp[5]=BBt[3]*B[2] + BBt[4]*B[5] + BBt[5]*B[8];

      tmp[6]=BBt[6]*B[0] + BBt[7]*B[3] + BBt[8]*B[6];
      tmp[7]=BBt[6]*B[1] + BBt[7]*B[4] + BBt[8]*B[7];
      tmp[8]=BBt[6]*B[2] + BBt[7]*B[5] + BBt[8]*B[8];

      // compute R as (a*B + 2*(l*adj(B)' - tmp))*denom; note that adj(B')=adj(B)'
      denom=l*(l*l-Bsq) - 2.0*detB;
      denom=1.0/denom;
      r(0)=(a*B[0] + 2.0*(l*adjB[0] - tmp[0]))*denom;
      r(1)=(a*B[1] + 2.0*(l*adjB[3] - tmp[1]))*denom;
      r(2)=(a*B[2] + 2.0*(l*adjB[6] - tmp[2]))*denom;

      r(3)=(a*B[3] + 2.0*(l*adjB[1] - tmp[3]))*denom;
      r(4)=(a*B[4] + 2.0*(l*adjB[4] - tmp[4]))*denom;
      r(5)=(a*B[5] + 2.0*(l*adjB[7] - tmp[5]))*denom;

      r(6)=(a*B[6] + 2.0*(l*adjB[2] - tmp[6]))*denom;
      r(7)=(a*B[7] + 2.0*(l*adjB[5] - tmp[7]))*denom;
      r(8)=(a*B[8] + 2.0*(l*adjB[8] - tmp[8]))*denom;
      }

      //double R[9];
      //r=Eigen::Map<Eigen::Matrix<double, 9, 1>>(R);
    }

    //
    // Produce a distance from being orthogonal for a random 3x3 matrix
    // Matrix is provided as a vector
    inline static double OrthogonalityError(const Eigen::Matrix<double, 9, 1>& a)
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
    // Compute the 3D null space (N) and 6D normal space (H) of the constraint Jacobian 
    // at a 9D vector r (not necessarilly a rotation-yet it should be rank-3)
    static void RowAndNullSpace(const Eigen::Matrix<double, 9, 1>& r, 
				      Eigen::Matrix<double, 9, 6>& H, 
				      Eigen::Matrix<double, 9, 3>& N,
				      Eigen::Matrix<double, 6, 6>& K,
				const double& norm_threshold = 0.1 ); 
   
  };
  

}

#endif
