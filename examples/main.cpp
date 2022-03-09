#include <iostream>
#include <cmath>
#include <sqpnp.h>
#include <opencv2/core.hpp>
#include <vector>
#include <chrono>
#include <assert.h>
#include <random>

//
// Generate noisy PnP data
void GenerateSyntheticPoints(int n, 
			     cv::Matx<double, 3, 3>& R,
			     cv::Vec<double, 3>& t,
			     std::vector<cv::Point3_<double>>& points,
			     std::vector<cv::Point_<double>>& projections,
			     std::vector<cv::Point_<double>>& noisy_projections,
			     const double& std_pixel_noise = 0.0,
			     const double& radius = 3.0 
			    )
{
    assert( n > 2);
    
    cv::Matx<double, 3,3> K(1400, 0, 1000, 
			    0, 1400, 900,
			    0, 0, 1);

    const double std_noise = std_pixel_noise / 1400;
    const double depth = 7*radius; // depth of the barycenter of the points
    
    const cv::Point3_<double> C(radius / 4, radius / 4, depth );
    
    
    // Generate a rotation matrix near the origin
    cv::Vec<double, 3> psi;// = mvnrnd([0; 0; 0], 0.001 * eye(3))';
    
    static std::random_device r;
    static std::default_random_engine generator(r());
    double sigma_psi = 0.1;
    std::normal_distribution<double> psi_noise(0.0, sigma_psi);
    psi[0] = psi_noise(generator);
    psi[1] = psi_noise(generator);
    psi[2] = psi_noise(generator);
    
    double sq_norm_psi = psi[0]*psi[0] + psi[1]*psi[1] + psi[2]*psi[2];
    double inv_w = 1.0 / (1 + sq_norm_psi );
    double s = (1 - sq_norm_psi) * inv_w,
	   v1 = 2*psi[0] * inv_w,
	   v2 = 2*psi[1] * inv_w,
	   v3 = 2*psi[2] * inv_w;
    R(0, 0)  = s*s + v1*v1 - v2*v2 - v3*v3;    R(0, 1) = 2*( v1*v2 - s*v3 );              R(0, 2) = 2*( v1*v3 + s*v2 );
    R(1, 0) = 2*( v1*v2 + s*v3 );              R(1, 1) = s*s - v1*v1 + v2*v2 - v3*v3;     R(1, 2) = 2*( v2*v3 - s*v1 ); 
    R(2, 0) = 2*( v1*v3 - s*v2 );              R(2, 1) = 2*( v2*v3 + s*v1 );              R(2, 2) = s*s - v1*v1 - v2*v2 + v3*v3;
    
    // Generate a translation that's about 1/25 of the depth
    std::normal_distribution<double> camera_position(0.0, depth/25 );
    cv::Vec<double, 3> pos(camera_position(generator), camera_position(generator), camera_position(generator));
    
    points.clear();
    projections.clear();
    noisy_projections.clear();
    std::normal_distribution<double> projection_noise(0.0, std_noise );
    
    while(points.size() < n)
    {
      std::normal_distribution<double> point_X(C.x,  radius); 
      std::normal_distribution<double> point_Y(C.y,  radius); 
      std::normal_distribution<double> point_Z(C.z,  radius); 
      
      cv::Vec<double, 3> Mw(point_X(generator), point_Y(generator), point_Z(generator)  );
      cv::Vec<double, 3> Mc = R*(Mw - pos);
       if ( Mc[2] < 0 )
       {
        std::cout << " Negative depth! Skipping...\n";
	continue;
       }
       cv::Vec<double, 2> proj( Mc[0] / Mc[2], Mc[1] / Mc[2] );
       // Add noise to projection
       cv::Vec<double, 2> noisy_proj = proj;
       noisy_proj[0] += projection_noise(generator);
       noisy_proj[1] += projection_noise(generator);
       noisy_proj[2] += projection_noise(generator);
       
       points.push_back(Mw);
       projections.push_back(proj);
       noisy_projections.push_back(noisy_proj);
    }
    
    t = -R*pos;
}

// compute the translational and angular error between the pose R,t contained in an SQPnP solution and the (true) pose Rg,tg
static void poseError(const sqpnp::SQPSolution& solution, const cv::Matx<double, 3, 3>& Rg, const cv::Vec<double, 3>& tg, double& terr, double& aerr)
{
    // translational error
    double a = tg(0) - solution.t(0);
    double b = tg(1) - solution.t(1);
    double c = tg(2) - solution.t(2);
    terr = sqrt(a*a + b*b + c*c);

    /* angular error, defined as the amount of rotation about a unit vector that transfers Rg to R.
     * The (residual) angle is computed with the inverse Rodrigues rotation formula
     */

    // compute trc as the trace of Rg'*R
    a = Rg(0, 0)*solution.r_hat[0] + Rg(1, 0)*solution.r_hat[3] + Rg(2, 0)*solution.r_hat[6];
    b = Rg(0, 1)*solution.r_hat[1] + Rg(1, 1)*solution.r_hat[4] + Rg(2, 1)*solution.r_hat[7];
    c = Rg(0, 2)*solution.r_hat[2] + Rg(1, 2)*solution.r_hat[5] + Rg(2, 2)*solution.r_hat[8];
    const double trc = a + b + c;
    a = 0.5*(trc - 1.0);
    aerr = acos(std::min(std::max(-1.0, a), 1.0)); // clamp to [-1, 1]
}


int main()
{
  int N = 105;
  int n = 10;
  double std_pixels = sqrt(7);
  
  std::vector<std::vector<cv::Point3_<double>>> vpoints;  
  std::vector<std::vector<cv::Point_<double>>> vprojections;
  std::vector<std::vector<cv::Point_<double>>> vnoisy_projections;
  
  std::vector<cv::Matx<double, 3, 3>> vRt;
  std::vector<cv::Vec<double, 3>> vtt;
  
  for (int i = 0; i < N; i++)
  {
    cv::Matx<double, 3, 3> Rt;
    cv::Vec<double, 3> tt;
    std::vector<cv::Point3_<double>> points;
    std::vector<cv::Point_<double>> projections;
    std::vector<cv::Point_<double>> noisy_projections;
    
    GenerateSyntheticPoints(n, Rt, tt, points, projections, noisy_projections, std_pixels);
    
    vpoints.push_back(points);
    vprojections.push_back(projections);
    vnoisy_projections.push_back(noisy_projections);
    vRt.push_back(Rt);
    vtt.push_back(tt);
  }
  
  auto start = std::chrono::steady_clock::now();
  double max_sq_error = 0, max_sq_proj_error = 0;
  std::vector<sqpnp::SQPSolution> solutions;
  for (int i = 0; i < N; i++)
  {
    sqpnp::PnPSolver solver(vpoints.at(i), vnoisy_projections.at(i));
    if (solver.IsValid() ) 
    {
      solver.Solve();
      if ( max_sq_error < solver.SolutionPtr(0)->sq_error ) 
      {
	   max_sq_error = solver.SolutionPtr(0)->sq_error;
           max_sq_proj_error = solver.AverageSquaredProjectionErrors().at(0);
      }
      solutions.push_back(*solver.SolutionPtr(0));
    }
  }
  
  auto finish = std::chrono::steady_clock::now();
  
  for (int i = 0; i < N; i++)
  {
    double terr, aerr;

    std::cout << i << "-th Solution : " << solutions.at(i);
    std::cout << i << "-th Rt : " << vRt.at(i) << std::endl;
    std::cout << i << "-th tt : " << vtt.at(i) << "\n";

    poseError(solutions.at(i), vRt.at(i), vtt.at(i), terr, aerr);
    std::cout << i << " translational error : " << terr <<  "  angular error : " << aerr*180.0/3.14159 << "\n\n";
  }
  
  auto diff = finish - start;
  std::cout << " Average execution time : " << std::chrono::duration_cast<std::chrono::microseconds> (diff).count() / N << std::endl;
  std::cout << " Maximum squared error : " << max_sq_error << std::endl;
  std::cout << " Maximum average squared projection error : " << max_sq_proj_error << std::endl;
  
  return 1;
}
