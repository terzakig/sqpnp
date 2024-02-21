//
// sqpnp.cpp
//
// George Terzakis (terzakig-at-hotmail-dot-com), September 2020
// Optimizations by Manolis Lourakis, February 2022, February 2024
// 
// Implementation of SQPnP as described in the paper:
//
// "A Consistently Fast and Globally Optimal Solution to the Perspective-n-Point Problem" by G. Terzakis and M. Lourakis
//     a) Paper:         https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460.pdf
//     b) Supplementary: https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460-supp.pdf

#include <sqpnp.h>

namespace sqpnp 
{
  
  const double SolverParameters::DEFAULT_RANK_TOLERANCE = 1e-7;
  const double SolverParameters::DEFAULT_SQP_SQUARED_TOLERANCE = 1e-10;
  const double SolverParameters::DEFAULT_SQP_DET_THRESHOLD = 1.001;
  const OmegaNullspaceMethod  SolverParameters::DEFAULT_OMEGA_NULLSPACE_METHOD = OmegaNullspaceMethod::RRQR;
  const NearestRotationMethod SolverParameters::DEFAULT_NEAREST_ROTATION_METHOD = NearestRotationMethod::FOAM;
  const double SolverParameters::DEFAULT_ORTHOGONALITY_SQUARED_ERROR_THRESHOLD = 1e-8;
  const double SolverParameters::DEFAULT_EQUAL_VECTORS_SQUARED_DIFF = 1e-10;
  const double SolverParameters::DEFAULT_EQUAL_SQUARED_ERRORS_DIFF = 1e-6;
  const double SolverParameters::DEFAULT_POINT_VARIANCE_THRESHOLD = 1e-5;
  
    
  const double PnPSolver::SQRT3 = std::sqrt(3.0);
  
  
  void PnPSolver::HandleSolution(SQPSolution& solution, double& min_sq_error)
  {
      bool cheirok = TestPositiveDepth( solution ) || TestPositiveMajorityDepths ( solution ); // check the majority if the check with centroid fails
      if ( cheirok )
      {
	
	solution.sq_error = ( Omega_ * solution.r_hat ).dot( solution.r_hat );
	if ( fabs( min_sq_error - solution.sq_error ) > parameters_.equal_squared_errors_diff )
	{
	  if (min_sq_error > solution.sq_error)
	  {
	    min_sq_error = solution.sq_error;
	    solutions_[0] = solution;
	    num_solutions_ = 1;
	  }
	}
	else // look for a solution that's almost equal to this
	{
	  bool found = false;
	  for (int i = 0; i < num_solutions_; i++)
	  {
	    if ( ( solutions_[i].r_hat - solution.r_hat ).squaredNorm() < parameters_.equal_vectors_squared_diff )
	    {
	      if ( solutions_[i].sq_error > solution.sq_error ) 
	      {
		solutions_[i] = solution;
	      }
	      found = true;
	      break;
	    }
	  }
	  if (!found)
	  {
	    solutions_[num_solutions_++] = solution;
	  }
	  if ( min_sq_error > solution.sq_error ) min_sq_error = solution.sq_error;
	}
      }
  }
  
  
  //
  // Solve the PnP 
  bool PnPSolver::Solve()
  {
    if ( !flag_valid_ ) return false;
  
    double min_sq_error = std::numeric_limits<double>::max();
    int num_eigen_points = num_null_vectors_ > 0 ? num_null_vectors_ : 1;
    // clear solutions
    num_solutions_ = 0;
    
    for (int i = 9 - num_eigen_points; i < 9; i++) 
    {
      // NOTE: No need to scale by sqrt(3) here, but better be there for other computations (i.e., orthogonality test)
      const Eigen::Matrix<double, 9, 1> e = SQRT3 * Eigen::Map<Eigen::Matrix<double, 9, 1>>( U_.block<9, 1>(0, i).data() );
      double orthogonality_sq_error = OrthogonalityError(e);
      // Find nearest rotation vector
      SQPSolution solution[2];
      
      // Avoid SQP if e is orthogonal
      if ( orthogonality_sq_error < parameters_.orthogonality_squared_error_threshold ) 
      {
	 solution[0].r_hat = Determinant9x1(e) * e;
	 solution[0].t = P_*solution[0].r_hat;
	 solution[0].num_iterations = 0;
	 
	 HandleSolution( solution[0], min_sq_error );
      }
      else
      {
	  NearestRotationMatrix( e, solution[0].r );
	  solution[0] = RunSQP( solution[0].r );
	  solution[0].t = P_*solution[0].r_hat;
	  HandleSolution( solution[0] , min_sq_error );

	  NearestRotationMatrix( -e, solution[1].r );
	  solution[1] = RunSQP( solution[1].r );
	  solution[1].t = P_*solution[1].r_hat;
	  HandleSolution( solution[1] , min_sq_error );
      }
    }

    int index, c = 1;
    while ((index = 9 - num_eigen_points - c) > 0 && min_sq_error > 3 * s_[index]) 
    {      
      const Eigen::Matrix<double, 9, 1> e = Eigen::Map<Eigen::Matrix<double, 9, 1>>( U_.block<9, 1>(0, index).data() );
      SQPSolution solution[2];
      
	NearestRotationMatrix( e, solution[0].r);
	solution[0] = RunSQP( solution[0].r );
	solution[0].t = P_*solution[0].r_hat;
	HandleSolution( solution[0], min_sq_error );

	NearestRotationMatrix( -e, solution[1].r);
	solution[1] = RunSQP( solution[1].r );
	solution[1].t = P_*solution[1].r_hat;
	HandleSolution( solution[1], min_sq_error );

      c++;
    }

    return true;
  }

  
  //
  // Run sequential quadratic programming on orthogonal matrices
  SQPSolution PnPSolver::RunSQP(const Eigen::Matrix<double, 9, 1>& r0)
  {
    Eigen::Matrix<double, 9, 1> r = r0;
    
    double delta_squared_norm = std::numeric_limits<double>::max();
    Eigen::Matrix<double, 9, 1> delta;
    int step = 0;
    
    while ( delta_squared_norm > parameters_.sqp_squared_tolerance && step++ < parameters_.sqp_max_iteration )
    {    
      SolveSQPSystem(r, delta);
      r += delta;
      delta_squared_norm = delta.squaredNorm();
    }
    
    SQPSolution solution;
    solution.num_iterations = step;
    solution.r = r;
    // clear the estimate and/or flip the matrix sign if necessary
    double det_r = Determinant9x1(solution.r);
    if (det_r < 0) 
    {
      solution.r = -r;
      det_r = -det_r;
    }
    if ( det_r > parameters_.sqp_det_threshold )
    {
      NearestRotationMatrix( solution.r, solution.r_hat );
    }
    else
    {
      solution.r_hat = solution.r;
    }
      
    
    return solution;
  }

  // Solve A*x=b for 3x3 SPD A.
  // The solution involves computing a lower triangular sqrt-free Cholesky factor
  // A=L*D*L' (L has ones on its diagonal, D is diagonal).
  //
  // Only the lower triangular part of A is accessed.
  //
  // The function returns 0 if successful, non-zero otherwise
  //
  // see http://euler.nmt.edu/~brian/ldlt.html
  //
  inline static int AxbSolveLDLt3x3(const Eigen::Matrix<double, 3, 3>& A, const Eigen::Matrix<double, 3, 1>& b,
                                    Eigen::Matrix<double, 3, 1>& x)
  {
    double L[3 * 3], v[2];

    // D is stored in L's diagonal, i.e. L[0], L[4], L[8]
    // its elements should be positive
    v[0] = L[0] = A(0, 0);
    if (v[0] < 1E-10) return 1;
    v[1] = 1.0 / v[0];
    L[3] = A(1, 0) * v[1];
    L[6] = A(2, 0) * v[1];
    // L[1] = L[2] = 0.0;

    v[0] = L[3] * L[0];
    v[1] = L[4] = A(1, 1) - L[3] * v[0];
    if (v[1] < 1E-10) return 2;
    L[7] = (A(2, 1) - L[6] * v[0]) / v[1];
    // L[5] = 0.0;

    v[0] = L[6] * L[0];
    v[1] = L[7] * L[4];
    L[8] = A(2, 2) - L[6] * v[0] - L[7] * v[1];

    // Forward solve L*x = b
    x[0] = b[0];
    x[1] = b[1] - L[3] * x[0];
    x[2] = b[2] - L[6] * x[0] - L[7] * x[1];

    // Backward solve D*L'*x = y
    x[2] = x[2] / L[8];
    x[1] = x[1] / L[4] - L[7] * x[2];
    x[0] = x[0] / L[0] - L[3] * x[1] - L[6] * x[2];

    return 0;
  }

  //
  // Solve the SQP system efficiently
  void PnPSolver::SolveSQPSystem(const Eigen::Matrix<double, 9, 1>& r, Eigen::Matrix<double, 9, 1>& delta)
  {
    double sqnorm_r1 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2], 
	   sqnorm_r2 = r[3]*r[3] + r[4]*r[4] + r[5]*r[5], 
	   sqnorm_r3 = r[6]*r[6] + r[7]*r[7] + r[8]*r[8];
    double dot_r1r2 = r[0]*r[3] + r[1]*r[4] + r[2]*r[5], 
	   dot_r1r3 = r[0]*r[6] + r[1]*r[7] + r[2]*r[8], 
	   dot_r2r3 = r[3]*r[6] + r[4]*r[7] + r[5]*r[8];
    
    // Obtain 6D normal (H) and 3D null space of the constraint Jacobian-J at the estimate (r)
    // NOTE: This is done via Gram-Schmidt orthogonalization
    Eigen::Matrix<double, 9, 3> N;  // Null space of J
    Eigen::Matrix<double, 9, 6> H;  // Row space of J
    Eigen::Matrix<double, 6, 6> JH; // The lower triangular matrix J*Q
    
    RowAndNullSpace(r, H, N, JH);
    
    // Great, now if delta = H*x + N*y, we first compute x by solving:
    // 
    //              (J*H)*x = g
    //
    // where g is the constraint vector g = [   1 - norm(r1)^2;
    // 					     	   1 - norm(r2)^2;
    //					     	   1 - norm(r3)^2;
    //					           -r1'*r2; 
    //						   -r2'*r3; 
    //						   -r1'*r3 ];
    Eigen::Matrix<double, 6, 1> g; 
    g[0] = 1 - sqnorm_r1; g[1] = 1 - sqnorm_r2; g[2] = 1 - sqnorm_r3; g[3] = -dot_r1r2; g[4] = -dot_r2r3; g[5] = -dot_r1r3;
      
    Eigen::Matrix<double, 6, 1> x;
    x[0] = g[0] / JH(0, 0);
    x[1] = g[1] / JH(1, 1);
    x[2] = g[2] / JH(2, 2);
    x[3] = ( g[3] - JH(3, 0)*x[0] - JH(3, 1)*x[1] ) / JH(3, 3);
    x[4] = ( g[4] - JH(4, 1)*x[1] - JH(4, 2)*x[2] - JH(4, 3)*x[3] ) / JH(4, 4);
    x[5] = ( g[5] - JH(5, 0)*x[0] - JH(5, 2)*x[2] - JH(5, 3)*x[3] - JH(5, 4)*x[4] ) / JH(5, 5);
    
    // Now obtain the component of delta in the row space of E as delta_h = H*x and assign straight into delta
    delta = H * x;
    
    // Then, solve for y from W*y = ksi, where matrix W and vector ksi are :
    //
    // W = N'*Omega*N and ksi = -N'*Omega*( r + delta_h );
    Eigen::Matrix<double, 3, 9> NtOmega = N.transpose() * Omega_ ;
    Eigen::Matrix<double, 3, 3> W = NtOmega * N;
    Eigen::Matrix<double, 3, 1> y, rhs = -NtOmega * ( delta + r );

    // solve with LDLt and if it fails, use LU
    if (AxbSolveLDLt3x3(W, rhs, y))
    {
      y = W.lu().solve(rhs);
    }

    // Finally, accumulate delta with component in tangent space (delta_n)
    delta += N*y;
  }
  
  
  //
  // Compute the 3D null space (N) and 6D normal space (H) of the constraint Jacobian at a 9D vector r 
  // (r is not necessarily a rotation but it must represent a rank-3 matrix)
  // NOTE: K is lower-triangular, so upper triangle may contain trash (is not filled by the function)...
  void PnPSolver::RowAndNullSpace(const Eigen::Matrix<double, 9, 1>& r, 
				    Eigen::Matrix<double, 9, 6>& H, // Row space 
				    Eigen::Matrix<double, 9, 3>& N, // Null space
				    Eigen::Matrix<double, 6, 6>& K,  // J*Q (J - Jacobian of constraints)
				  const double& norm_threshold // Used to discard columns of Pn when finding null space
 				) // threshold for column vector norm (of Pn)
  {
    // Applying Gram-Schmidt orthogonalization on the Jacobian. 
    // The steps are fixed here to take advantage of the sparse form of the matrix
    //
    H = Eigen::Matrix<double, 9, 6>::Zero();
    
    // 1. q1
    double norm_r1 = sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] );
    double inv_norm_r1 = norm_r1 > 1e-5 ? 1.0 / norm_r1 : 0.0;
    H(0, 0) = r[0] * inv_norm_r1; H(1, 0) = r[1] * inv_norm_r1; H(2, 0) = r[2] * inv_norm_r1;
    K(0, 0) = 2*norm_r1;
    
    // 2. q2 
    double norm_r2 = sqrt( r[3]*r[3] + r[4]*r[4] + r[5]*r[5] );
    double inv_norm_r2 = 1.0 / norm_r2;
    H(3, 1) = r[3]*inv_norm_r2; H(4, 1) = r[4]*inv_norm_r2; H(5, 1) = r[5]*inv_norm_r2;
    K(1, 0) = 0; K(1, 1) = 2*norm_r2;
    
    // 3. q3 = (r3'*q2)*q2 - (r3'*q1)*q1 ; q3 = q3/norm(q3)
    double norm_r3 = sqrt( r[6]*r[6] + r[7]*r[7] + r[8]*r[8] );
    double inv_norm_r3 = 1.0 / norm_r3;
    H(6, 2) = r[6]*inv_norm_r3; H(7, 2) = r[7]*inv_norm_r3; H(8, 2) = r[8]*inv_norm_r3;
    K(2, 0) = K(2, 1) = 0; K(2, 2) = 2*norm_r3;
    
    // 4. q4
    double dot_j4q1 = r[3]*H(0, 0) + r[4]*H(1, 0) + r[5]*H(2, 0),
	   dot_j4q2 = r[0]*H(3, 1) + r[1]*H(4, 1) + r[2]*H(5, 1);
    
    H(0, 3) = r[3] - dot_j4q1*H(0, 0); H(1, 3) = r[4] - dot_j4q1*H(1, 0); H(2, 3) = r[5] - dot_j4q1*H(2, 0);
    H(3, 3) = r[0] - dot_j4q2*H(3, 1); H(4, 3) = r[1] - dot_j4q2*H(4, 1); H(5, 3) = r[2] - dot_j4q2*H(5, 1);
    double inv_norm_j4 = 1.0 / sqrt( H(0, 3)*H(0, 3) + H(1, 3)*H(1, 3) + H(2, 3)*H(2, 3) + 
				     H(3, 3)*H(3, 3) + H(4, 3)*H(4, 3) + H(5, 3)*H(5, 3) );
    
    H(0, 3) *= inv_norm_j4; H(1, 3) *= inv_norm_j4; H(2, 3) *= inv_norm_j4;
    H(3, 3) *= inv_norm_j4; H(4, 3) *= inv_norm_j4; H(5, 3) *= inv_norm_j4;
    
    K(3, 0) = r[3]*H(0, 0) + r[4]*H(1, 0) + r[5]*H(2, 0); K(3, 1) = r[0]*H(3, 1) + r[1]*H(4, 1) + r[2]*H(5, 1); 
    K(3, 2) = 0; K(3, 3) = r[3]*H(0, 3) + r[4]*H(1, 3) + r[5]*H(2, 3)  +  r[0]*H(3, 3) + r[1]*H(4, 3) + r[2]*H(5, 3);
    
    // 5. q5
    double dot_j5q2 = r[6]*H(3, 1) + r[7]*H(4, 1) + r[8]*H(5, 1),
	   dot_j5q3 = r[3]*H(6, 2) + r[4]*H(7, 2) + r[5]*H(8, 2),
	   dot_j5q4 = r[6]*H(3, 3) + r[7]*H(4, 3) + r[8]*H(5, 3);
    
    H(0, 4) = -dot_j5q4*H(0, 3); 			  H(1, 4) = -dot_j5q4*H(1, 3); 			H(2, 4) = -dot_j5q4*H(2, 3);
    H(3, 4) = r[6] - dot_j5q2*H(3, 1) - dot_j5q4*H(3, 3); H(4, 4) = r[7] - dot_j5q2*H(4, 1) - dot_j5q4*H(4, 3); H(5, 4) = r[8] - dot_j5q2*H(5, 1) - dot_j5q4*H(5, 3);
    H(6, 4) = r[3] - dot_j5q3*H(6, 2); H(7, 4) = r[4] - dot_j5q3*H(7, 2); H(8, 4) = r[5] - dot_j5q3*H(8, 2);
    
    H.block<9, 1>(0, 4) *= (1.0 / H.col(4).norm());
   
    K(4, 0) = 0; K(4, 1) = r[6]*H(3, 1) + r[7]*H(4, 1) + r[8]*H(5, 1); K(4, 2) = r[3]*H(6, 2) + r[4]*H(7, 2) + r[5]*H(8, 2);
    K(4, 3) = r[6]*H(3, 3) + r[7]*H(4, 3) + r[8]*H(5, 3); 
    K(4, 4) = r[6]*H(3, 4) + r[7]*H(4, 4) + r[8]*H(5, 4)  +  r[3]*H(6, 4) + r[4]*H(7, 4) + r[5]*H(8, 4); 
    
    // 4. q6
    double dot_j6q1 = r[6]*H(0, 0) + r[7]*H(1, 0) + r[8]*H(2, 0),
	   dot_j6q3 = r[0]*H(6, 2) + r[1]*H(7, 2) + r[2]*H(8, 2), 
	   dot_j6q4 = r[6]*H(0, 3) + r[7]*H(1, 3) + r[8]*H(2, 3), 
	   dot_j6q5 = r[0]*H(6, 4) + r[1]*H(7, 4) + r[2]*H(8, 4)  +  r[6]*H(0, 4) + r[7]*H(1, 4) + r[8]*H(2, 4);
    
    H(0, 5) = r[6] - dot_j6q1*H(0, 0) - dot_j6q4*H(0, 3) - dot_j6q5*H(0, 4); 
    H(1, 5) = r[7] - dot_j6q1*H(1, 0) - dot_j6q4*H(1, 3) - dot_j6q5*H(1, 4); 
    H(2, 5) = r[8] - dot_j6q1*H(2, 0) - dot_j6q4*H(2, 3) - dot_j6q5*H(2, 4);
    
    H(3, 5) = -dot_j6q5*H(3, 4) - dot_j6q4*H(3, 3); 
    H(4, 5) = -dot_j6q5*H(4, 4) - dot_j6q4*H(4, 3); 
    H(5, 5) = -dot_j6q5*H(5, 4) - dot_j6q4*H(5, 3);
    
    H(6, 5) = r[0] - dot_j6q3*H(6, 2) - dot_j6q5*H(6, 4); 
    H(7, 5) = r[1] - dot_j6q3*H(7, 2) - dot_j6q5*H(7, 4); 
    H(8, 5) = r[2] - dot_j6q3*H(8, 2) - dot_j6q5*H(8, 4);
    
    H.block<9, 1>(0, 5) *= (1.0 / H.col(5).norm());
    
    K(5, 0) = r[6]*H(0, 0) + r[7]*H(1, 0) + r[8]*H(2, 0); K(5, 1) = 0; K(5, 2) = r[0]*H(6, 2) + r[1]*H(7, 2) + r[2]*H(8, 2);
    K(5, 3) = r[6]*H(0, 3) + r[7]*H(1, 3) + r[8]*H(2, 3); K(5, 4) = r[6]*H(0, 4) + r[7]*H(1, 4) + r[8]*H(2, 4) +   r[0]*H(6, 4) + r[1]*H(7, 4) + r[2]*H(8, 4);
    K(5, 5) = r[6]*H(0, 5) + r[7]*H(1, 5) + r[8]*H(2, 5) + r[0]*H(6, 5) + r[1]*H(7, 5) + r[2]*H(8, 5);
    
    // Great! Now H is an orthogonalized, sparse basis of the Jacobian row space and K is filled.
    //
    // Now get a projector onto the null space of H:
    const Eigen::Matrix<double, 9, 9> Pn = Eigen::Matrix<double, 9, 9>::Identity() - ( H*H.transpose() );
    
    // Now we need to pick 3 columns of P with non-zero norm (> 0.3) and some angle between them (> 0.3).
    //
    // Find the 3 columns of Pn with largest norms
    int index1 = -1, 
	index2 = -1, 
	index3 = -1;
    double  max_norm1 = std::numeric_limits<double>::min(),
	    min_dot12 = std::numeric_limits<double>::max(),
	    min_dot1323 = std::numeric_limits<double>::max();
    
    
    double col_norms[9];
    for (int i = 0; i < 9; i++)
    {
      col_norms[i] = Pn.col(i).norm();
      if ( col_norms[i] >= norm_threshold)
      {
	if (max_norm1 < col_norms[i])
	{
	  max_norm1 = col_norms[i];
	  index1 = i;
	}
      }
    }
    const auto& v1 = Pn.block<9, 1>(0, index1);
    N.block<9, 1>(0, 0) = v1 * ( 1.0 / max_norm1 );
    col_norms[index1] = -1.0; // mark to avoid use in subsequent loops
    
    for (int i = 0; i < 9; i++)
    {
      //if (i == index1) continue;
      if ( col_norms[i] >= norm_threshold)
      {
	double cos_v1_x_col = fabs(Pn.col(i).dot(v1) / col_norms[i]);
	
	if ( cos_v1_x_col <= min_dot12)
	{
	  index2 = i;
	  min_dot12 = cos_v1_x_col;
	}
      }
    }
    const auto& v2 = Pn.block<9, 1>(0, index2);
    N.block<9, 1>(0, 1) = v2 - v2.dot( N.col(0) ) * N.col(0);
    N.block<9, 1>(0, 1) *= (1.0 / N.col(1).norm());
    col_norms[index2] = -1.0; // mark to avoid use in loop below
    
    for (int i = 0; i < 9; i++)
    {
      //if (i == index2 || i == index1) continue;
      if ( col_norms[i] >= norm_threshold)
      {
	double inv_norm = 1.0 / col_norms[i];
	double cos_v1_x_col = fabs(Pn.col(i).dot(v1) * inv_norm);
	double cos_v2_x_col = fabs(Pn.col(i).dot(v2) * inv_norm);
	
	if ( cos_v1_x_col + cos_v2_x_col <= min_dot1323)
	{
	  index3 = i;
	  min_dot1323 = cos_v2_x_col + cos_v2_x_col;
	}
      }
    }
    
    // Now orthogonalize the remaining 2 vectors v2, v3 into N
    const auto& v3 = Pn.block<9, 1>(0, index3);
    
    N.block<9, 1>(0, 2) = v3 - ( v3.dot( N.col(1) ) * N.col(1) ) - ( v3.dot( N.col(0) ) * N.col(0) );
    N.block<9, 1>(0, 2) *= (1.0 / N.col(2).norm());
    
  }
  
}
