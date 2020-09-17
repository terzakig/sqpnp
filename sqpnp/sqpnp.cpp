//
// sqpnp.cpp
//
// George Terzakis (terzakig-at-hotmail-dot-com), September 2020
// 
// Implementation of SQPnP as described in the paper:
//
// "A Consistently Fast and Globally Optimal Solution to the Perspective-n-Point Problem" by G. Terzakis and M. Lourakis
//  	 a) Paper: 	   http://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460.pdf 
//       b) Supplementary: https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460.pdf

#include <sqpnp.h>

namespace sqpnp 
{
  
  const double SolverParameters::DEFAULT_RANK_TOLERANCE = 1e-7;
  const double SolverParameters::DEFAULT_SQP_SQUARED_TOLERANCE = 1e-10;
  const double SolverParameters::DEFAULT_SQP_DET_THRESHOLD = 1.001;
  const double SolverParameters::DEFAULT_ORTHOGONALITY_SQUARED_ERROR_THRESHOLD = 1e-8;
  const double SolverParameters::DEFAULT_EQUAL_VECTORS_SQUARED_DIFF = 1e-10;
  const double SolverParameters::DEFAULT_EQUAL_SQUARED_ERRORS_DIFF = 1e-6;
    
  const double PnPSolver::SQRT3 = std::sqrt(3);
  
  
  void PnPSolver::HandleSolution(SQPSolution& solution, double& min_sq_error)
  {
      if ( TestPositiveDepth( solution ) )
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
      const Eigen::Vector<double, 9> e = SQRT3 * U_.block<9, 1>(0, i).reshaped(9, 1);
      double orthogonality_sq_error = OrthogonalityError(e);
      // Find nearest rotation vector
      SQPSolution solution[2];
      
      // Avoid SQP if e is orthogonal
      if ( orthogonality_sq_error < parameters_.orthogonality_squared_error_threshold ) 
      {
	 solution[0].r_hat = Determinant3x3(e) * e;
	 solution[0].t = P_*solution[0].r_hat;
	 solution[0].num_iterations = 0;
	 
	 HandleSolution( solution[0], min_sq_error );
      }
      else
      {
	for (int k = 0; k < 2; k++)
	{
	  NearestRotationMatrix( (k == 0 ? 1 : -1) * e, solution[k].r );
	  solution[k] = RunSQP( solution[k].r );
	  solution[k].t = P_*solution[k].r_hat;
	  
	  HandleSolution( solution[k] , min_sq_error );
	}
      }
    }

    int c = 1;
    while (min_sq_error > 3 * s_[9 - num_eigen_points - c] && 9 - num_eigen_points - c > 0) 
    {      
      int index = 9 - num_eigen_points - c;

      const Eigen::Vector<double, 9> e = Eigen::Map<Eigen::Vector<double, 9>>( U_.block<9, 1>(0, index).data() );
      SQPSolution solution[2];
      
      for (int k = 0; k < 2; k++)
      {
	NearestRotationMatrix( (k == 0 ? 1 : -1)*e, solution[k].r);
      
	solution[k] = RunSQP( solution[k].r );
	solution[k].t = P_*solution[k].r_hat;
	
	HandleSolution( solution[k], min_sq_error );
      }
	
      c++;
    }

    return true;
  }

  
  //
  // Run sequential quadratic programming on orthogonal matrices
  SQPSolution PnPSolver::RunSQP(const Eigen::Vector<double, 9>& r0)
  {
    Eigen::Vector<double, 9> r = r0;
    
    double delta_squared_norm = std::numeric_limits<double>::max();
    Eigen::Vector<double, 9> delta;
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
    double det_r = Determinant3x3(solution.r);
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
  
  //
  // Solve the SQP system efficiently
  void PnPSolver::SolveSQPSystem(const Eigen::Vector<double, 9>& r, Eigen::Vector<double, 9>& delta )
  {
    double sqnorm_r1 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2], 
	   sqnorm_r2 = r[3]*r[3] + r[4]*r[4] + r[5]*r[5], 
	   sqnorm_r3 = r[6]*r[6] + r[7]*r[7] + r[8]*r[8];
    double dot_r1r2 = r[0]*r[3] + r[1]*r[4] + r[2]*r[5], 
	   dot_r1r3 = r[0]*r[6] + r[1]*r[7] + r[2]*r[8], 
	   dot_r2r3 = r[3]*r[6] + r[4]*r[7] + r[5]*r[8];
    
    // Obtain 6D normal (Q) and 3D null space of the constraint Jacobian-J at the estimate (r)
    // NOTE: Thsi is done via Gram-Schmidt orthogoalization
    Eigen::Matrix<double, 9, 3> N;  // Null space of J
    Eigen::Matrix<double, 9, 6> Q;  // Row space of J
    Eigen::Matrix<double, 6, 6> JQ; // The lower triangular matrix J*Q
    
    RowAndNullSpace(r, Q, N, JQ);
    
    // Great, now if delta = Q*x + N*y, we first compute x by solving:
    // 
    //              (J*Q)*x = g
    //
    // where g is the constraint vector g = [   1 - norm(r1)^2;
    // 					     	   1 - norm(r2)^2;
    //					     	   1 - norm(r3)^2;
    //					           -r1'*r2; 
    //						   -r2'*r3; 
    //						   -r1'*r3 ];
    Eigen::Vector<double, 6> g; 
    g[0] = 1 - sqnorm_r1; g[1] = 1 - sqnorm_r2; g[2] = 1 - sqnorm_r3; g[3] = -dot_r1r2; g[4] = -dot_r2r3; g[5] = -dot_r1r3;
      
    Eigen::Vector<double, 6> x;
    x[0] = g[0] / JQ(0, 0);
    x[1] = g[1] / JQ(1, 1);
    x[2] = g[2] / JQ(2, 2);
    x[3] = ( g[3] - JQ(3, 0)*x[0] - JQ(3, 1)*x[1] ) / JQ(3, 3);
    x[4] = ( g[4] - JQ(4, 1)*x[1] - JQ(4, 2)*x[2] - JQ(4, 3)*x[3] ) / JQ(4, 4);
    x[5] = ( g[5] - JQ(5, 0)*x[0] - JQ(5, 2)*x[2] - JQ(5, 3)*x[3] - JQ(5, 4)*x[4] ) / JQ(5, 5);
    
    // Now obtain the component of delta in the row space of E as delta_h = Q'*x and assign straint into delta
    delta = Q * x;
    
    // Finally, solve for y from W*y = ksi , where matrix W and vector ksi are :
    //
    // W = N'*Omega*N and ksi = -N'*Omega*( r + delta_h );
    Eigen::Matrix<double, 3, 9> NtOmega = N.transpose() * Omega_ ;
    Eigen::Matrix<double, 3, 3> W = NtOmega * N, Winv;
    InvertSymmetric3x3(W, Winv); // NOTE: This maybe also analytical with Eigen, but hey...
    
    Eigen::Vector<double, 3> y = -Winv * NtOmega * ( delta + r );
    
    // FINALLY, accumulate delta with component in tangent space (delta_n)
    delta += N*y;
  }
  
  
  //
  // Compute the 3D null space and 6D normal space of the constraint Jacobian at a 9D vector r 
  // (r is not necessarilly a rotation but it must represent an rank-3 matrix )
  // NOTE: K is lower-triangular, so upper triangle may contain trash (is not filled by the function)...
  void PnPSolver::RowAndNullSpace(const Eigen::Vector<double, 9>& r, 
				    Eigen::Matrix<double, 9, 6>& Q, // Row space 
				    Eigen::Matrix<double, 9, 3>& N, // Null space
				    Eigen::Matrix<double, 6, 6>& K,  // J*Q (J - Jacobian of constraints)
				  const double& norm_threshold // Used to discard columns of Pn when finding null space
 				) // threshold for column vector norm (of Pn)
  {
    // Applying Gram-Schmidt orthogonalization on the Jacobian. 
    // The steps are fixed here to take advantage of the sparse form of the matrix
    //
    Q = Eigen::Matrix<double, 9, 6>::Zero();
    
    // 1. q1
    double norm_r1 = sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] );
    double inv_norm_r1 = norm_r1 > 1e-5 ? 1.0 / norm_r1 : 0.0;
    Q(0, 0) = r[0] * inv_norm_r1; Q(1, 0) = r[1] * inv_norm_r1; Q(2, 0) = r[2] * inv_norm_r1;
    K(0, 0) = 2*norm_r1;
    
    // 2. q2 
    double norm_r2 = sqrt( r[3]*r[3] + r[4]*r[4] + r[5]*r[5] );
    double inv_norm_r2 = 1.0 / norm_r2;
    Q(3, 1) = r[3]*inv_norm_r2; Q(4, 1) = r[4]*inv_norm_r2; Q(5, 1) = r[5]*inv_norm_r2;
    K(1, 0) = 0; K(1, 1) = 2*norm_r2;
    
    // 3. q3 = (r3'*q2)*q2 - (r3'*q1)*q1 ; q3 = q3/norm(q3)
    double norm_r3 = sqrt( r[6]*r[6] + r[7]*r[7] + r[8]*r[8] );
    double inv_norm_r3 = 1.0 / norm_r3;
    Q(6, 2) = r[6]*inv_norm_r3; Q(7, 2) = r[7]*inv_norm_r3; Q(8, 2) = r[8]*inv_norm_r3;
    K(2, 0) = K(2, 1) = 0; K(2, 2) = 2*norm_r3;
    
    // 4. q4
    double dot_j4q1 = r[3]*Q(0, 0) + r[4]*Q(1, 0) + r[5]*Q(2, 0),
	   dot_j4q2 = r[0]*Q(3, 1) + r[1]*Q(4, 1) + r[2]*Q(5, 1);
    
    Q(0, 3) = r[3] - dot_j4q1*Q(0, 0); Q(1, 3) = r[4] - dot_j4q1*Q(1, 0); Q(2, 3) = r[5] - dot_j4q1*Q(2, 0);
    Q(3, 3) = r[0] - dot_j4q2*Q(3, 1); Q(4, 3) = r[1] - dot_j4q2*Q(4, 1); Q(5, 3) = r[2] - dot_j4q2*Q(5, 1);
    double inv_norm_j4 = 1.0 / sqrt( Q(0, 3)*Q(0, 3) + Q(1, 3)*Q(1, 3) + Q(2, 3)*Q(2, 3) + 
				     Q(3, 3)*Q(3, 3) + Q(4, 3)*Q(4, 3) + Q(5, 3)*Q(5, 3) );
    
    Q(0, 3) *= inv_norm_j4; Q(1, 3) *= inv_norm_j4; Q(2, 3) *= inv_norm_j4;
    Q(3, 3) *= inv_norm_j4; Q(4, 3) *= inv_norm_j4; Q(5, 3) *= inv_norm_j4;
    
    K(3, 0) = r[3]*Q(0, 0) + r[4]*Q(1, 0) + r[5]*Q(2, 0); K(3, 1) = r[0]*Q(3, 1) + r[1]*Q(4, 1) + r[2]*Q(5, 1); 
    K(3, 2) = 0; K(3, 3) = r[3]*Q(0, 3) + r[4]*Q(1, 3) + r[5]*Q(2, 3)  +  r[0]*Q(3, 3) + r[1]*Q(4, 3) + r[2]*Q(5, 3);
    
    // 5. q5
    double dot_j5q2 = r[6]*Q(3, 1) + r[7]*Q(4, 1) + r[8]*Q(5, 1),
	   dot_j5q3 = r[3]*Q(6, 2) + r[4]*Q(7, 2) + r[5]*Q(8, 2),
	   dot_j5q4 = r[6]*Q(3, 3) + r[7]*Q(4, 3) + r[8]*Q(5, 3);
    
    Q(0, 4) = -dot_j5q4*Q(0, 3); 			  Q(1, 4) = -dot_j5q4*Q(1, 3); 			Q(2, 4) = -dot_j5q4*Q(2, 3);
    Q(3, 4) = r[6] - dot_j5q2*Q(3, 1) - dot_j5q4*Q(3, 3); Q(4, 4) = r[7] - dot_j5q2*Q(4, 1) - dot_j5q4*Q(4, 3); Q(5, 4) = r[8] - dot_j5q2*Q(5, 1) - dot_j5q4*Q(5, 3);
    Q(6, 4) = r[3] - dot_j5q3*Q(6, 2); Q(7, 4) = r[4] - dot_j5q3*Q(7, 2); Q(8, 4) = r[5] - dot_j5q3*Q(8, 2);
    
    Q.block<9, 1>(0, 4) /= Q.col(4).norm();
   
    K(4, 0) = 0; K(4, 1) = r[6]*Q(3, 1) + r[7]*Q(4, 1) + r[8]*Q(5, 1); K(4, 2) = r[3]*Q(6, 2) + r[4]*Q(7, 2) + r[5]*Q(8, 2);
    K(4, 3) = r[6]*Q(3, 3) + r[7]*Q(4, 3) + r[8]*Q(5, 3); 
    K(4, 4) = r[6]*Q(3, 4) + r[7]*Q(4, 4) + r[8]*Q(5, 4)  +  r[3]*Q(6, 4) + r[4]*Q(7, 4) + r[5]*Q(8, 4); 
    
    
    // 4. q6
    double dot_j6q1 = r[6]*Q(0, 0) + r[7]*Q(1, 0) + r[8]*Q(2, 0),
	   dot_j6q3 = r[0]*Q(6, 2) + r[1]*Q(7, 2) + r[2]*Q(8, 2), 
	   dot_j6q4 = r[6]*Q(0, 3) + r[7]*Q(1, 3) + r[8]*Q(2, 3), 
	   dot_j6q5 = r[0]*Q(6, 4) + r[1]*Q(7, 4) + r[2]*Q(8, 4)  +  r[6]*Q(0, 4) + r[7]*Q(1, 4) + r[8]*Q(2, 4);
    
    Q(0, 5) = r[6] - dot_j6q1*Q(0, 0) - dot_j6q4*Q(0, 3) - dot_j6q5*Q(0, 4); 
    Q(1, 5) = r[7] - dot_j6q1*Q(1, 0) - dot_j6q4*Q(1, 3) - dot_j6q5*Q(1, 4); 
    Q(2, 5) = r[8] - dot_j6q1*Q(2, 0) - dot_j6q4*Q(2, 3) - dot_j6q5*Q(2, 4);
    
    Q(3, 5) = -dot_j6q5*Q(3, 4) - dot_j6q4*Q(3, 3); 
    Q(4, 5) = -dot_j6q5*Q(4, 4) - dot_j6q4*Q(4, 3); 
    Q(5, 5) = -dot_j6q5*Q(5, 4) - dot_j6q4*Q(5, 3);
    
    Q(6, 5) = r[0] - dot_j6q3*Q(6, 2) - dot_j6q5*Q(6, 4); 
    Q(7, 5) = r[1] - dot_j6q3*Q(7, 2) - dot_j6q5*Q(7, 4); 
    Q(8, 5) = r[2] - dot_j6q3*Q(8, 2) - dot_j6q5*Q(8, 4);
    
    Q.block<9, 1>(0, 5) /= Q.col(5).norm();
    
    K(5, 0) = r[6]*Q(0, 0) + r[7]*Q(1, 0) + r[8]*Q(2, 0); K(5, 1) = 0; K(5, 2) = r[0]*Q(6, 2) + r[1]*Q(7, 2) + r[2]*Q(8, 2);
    K(5, 3) = r[6]*Q(0, 3) + r[7]*Q(1, 3) + r[8]*Q(2, 3); K(5, 4) = r[6]*Q(0, 4) + r[7]*Q(1, 4) + r[8]*Q(2, 4) +   r[0]*Q(6, 4) + r[1]*Q(7, 4) + r[2]*Q(8, 4);
    K(5, 5) = r[6]*Q(0, 5) + r[7]*Q(1, 5) + r[8]*Q(2, 5) + r[0]*Q(6, 5) + r[1]*Q(7, 5) + r[2]*Q(8, 5);
    
    // Great! Now Q is an orthogonalized, sparse basis of the Jacobian row space and K is filled.
    //
    // Now get a projector onto the null space Q:
    const Eigen::Matrix<double, 9, 9> Pn = Eigen::Matrix<double,9, 9>::Identity() - ( Q*Q.transpose() ); 
    
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
    
    for (int i = 0; i < 9; i++)
    {
      if (i == index1) continue;
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
    N.block<9, 1>(0, 1) /= N.col(1).norm();
    
    for (int i = 0; i < 9; i++)
    {
      if (i == index2 || i == index1) continue;
      if ( col_norms[i] >= norm_threshold)
      {
	double cos_v1_x_col = fabs(Pn.col(i).dot(v1) / col_norms[i]);
	double cos_v2_x_col = fabs(Pn.col(i).dot(v2) / col_norms[i]);
	
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
    N.block<9, 1>(0, 2) /= N.col(2).norm();
    
    //std::cout << "========= indexes: " << index1 << " , " << index2 << " , " << index3 << std::endl;
  }
  
}