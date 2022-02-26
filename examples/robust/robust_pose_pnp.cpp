//
// robust_pose_pnp.cpp
//
// Manolis Lourakis (lourakis **at** ics forth gr), February 2022

#include <cmath>
#include <cstdio>
#include <iostream>
#include <cstring>
#include <vector>

#include <Eigen/Core>
#include <Eigen/Dense>

#include <RansacLib/ransac.h>
#include "robust_pose_pnp.h"

#undef ROBUST_POSE_PNP_DEBUG  // enables some sample size checks

namespace robust_pose_pnp {

PoseEstimator::PoseEstimator(const std::vector<sqpnp::_Point> *_3dpoints, const std::vector<sqpnp::_Projection> *_projections,
                             int samplesz, int nm_samplesz) {
  _3dpoints_ = _3dpoints;
  _projections_ = _projections;

  set_sample_sizes(samplesz, nm_samplesz);
}

// Robustly solves for camera pose lambda*x = R*X+t  with lambda > 0.
// Inlier and outlier indices are optionally returned (in increasing order) in idxInliers & idxOutliers.
// Returns 0 if successful, zero otherwise
int PoseEstimator::ransacfit(double inlPcent, double outlThresh, robust_pose_pnp::Matrix34d& best_pose, std::vector<int> *idxInliers, std::vector<int> *idxOutliers) const {
  ransac_lib::LORansacOptions options;

  options.min_num_iterations_ = 100u;
  options.max_num_iterations_ = 10000u;
  // make sure that min_num_iterations_ are at least as many as predicted for the given outlier ratio and sample size
  options.min_num_iterations_ = ransac_lib::utils::NumRequiredIterations(inlPcent, 1.0-0.99, samplesz_, options.min_num_iterations_, options.max_num_iterations_);
  options.squared_inlier_threshold_ = outlThresh * outlThresh;
  options.final_least_squares_ = true;

  // additional params with their defaults
  //options.min_sample_multiplicator_ = 7;
  //options.num_lsq_iterations_ = 4;
  //options.num_lo_steps_ = 10;

#if 0
  std::random_device rand_dev;
  options.random_seed_ = rand_dev();
#endif

  ransac_lib::LocallyOptimizedMSAC<robust_pose_pnp::Matrix34d,
                                  std::vector<robust_pose_pnp::Matrix34d, Eigen::aligned_allocator<robust_pose_pnp::Matrix34d> >,
                                  robust_pose_pnp::PoseEstimator>
  lomsac;
  ransac_lib::RansacStatistics ransac_stats;

//std::cout << "... running LOMSAC" << std::endl;
  int num_inliers = lomsac.EstimateModel(options, *this, &best_pose, &ransac_stats);

  int npts = _3dpoints_->size();
  int num_outliers = npts - num_inliers;

  if (idxInliers != nullptr) {
    idxInliers->clear();
    idxInliers->resize(num_inliers);

    for (int i=0; i<num_inliers; ++i)
      (*idxInliers)[i] = ransac_stats.inlier_indices[i];
  }

  if (idxOutliers != nullptr) {
    idxOutliers->clear();
    idxOutliers->resize(num_outliers);

    // num_inliers == ransac_stats.inlier_indices.size()
    std::vector<int> isoutl(npts, 1);
    for (int i=0; i<num_inliers; ++i)
      isoutl[ransac_stats.inlier_indices[i]] = 0;
    for (int i=0, j=0; i<npts; ++i)
      if (isoutl[i]) (*idxOutliers)[j++] = i;
  }

/****/
  std::cout << "... LOMSAC found " << num_inliers << " inliers in "
            << ransac_stats.num_iterations << " iterations with an inlier "
            << "ratio of " << ransac_stats.inlier_ratio << std::endl;
/****/

  return num_inliers == 0;
}

#if 0
// this version creates local copies for the points corresponding to sample
int PoseEstimator::sqpnp_solve(const std::vector<int>& sample,
                                 std::vector<robust_pose_pnp::Matrix34d, Eigen::aligned_allocator<robust_pose_pnp::Matrix34d> > *poses) const {
  const int nsample = sample.size();
//std::cout<< "in solve "<<nsample<<"\n"<<std::flush;

  // pose from sample
  std::vector<sqpnp::_Point> points(nsample);
  std::vector<sqpnp::_Projection> projections(nsample);

  for(int i=0; i<nsample; ++i){
    int j = sample[i];

    points[i] = (*_3dpoints_)[j];
    projections[i] = (*_projections_)[j];
  }

  sqpnp::SolverParameters params;
  // more relaxed convergence criteria below
  //params.sqp_max_iteration = 10;
  //params.rank_tolerance = 1E-06;
  //params.sqp_squared_tolerance = 1E-09;

  std::vector<double>wghts(nsample, 1.0);
  sqpnp::PnPSolver solver(points, projections, wghts, params);

  if(solver.IsValid()){
    solver.Solve();

    poses->resize(solver.NumberOfSolutions());
    for (int i = 0; i < solver.NumberOfSolutions(); i++)
    {
      const sqpnp::SQPSolution *sol = solver.SolutionPtr(i);

//      std::cout << "\nSolution " << i << ":\n";
//      std::cout << *sol << std::endl;
//      std::cout << " Average squared projection error : " << solver.AverageSquaredProjectionErrors().at(i) << std::endl;

      (*poses)[i].block<3,3>(0,0) = Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor> >( sol->r_hat.data() );
      (*poses)[i].block<3,1>(0,3) = Eigen::Map<const Eigen::Matrix<double, 3, 1, Eigen::ColMajor> >( sol->t.data() );
    }
  }
  else return 0;

  return solver.NumberOfSolutions();
}

#else
// this uses directly the correspondences using the weights: 1 for those in sample, 0 for the rest
int PoseEstimator::sqpnp_solve(const std::vector<int>& sample,
                                 std::vector<robust_pose_pnp::Matrix34d, Eigen::aligned_allocator<robust_pose_pnp::Matrix34d> > *poses) const {
  const int npts = _3dpoints_->size();
  const int nsample = sample.size();
//std::cout<< "in solve "<<nsample<<"\n"<<std::flush;

  // pose from sample
  std::vector<double>wghts(npts, 0.0);
  for(int i=0; i<nsample; ++i){
    wghts[sample[i]] = 1.0;
  }

  sqpnp::SolverParameters params;
  // more relaxed convergence parameters next
  //params.sqp_max_iteration = 10;
  //params.rank_tolerance = 1E-06;
  //params.sqp_squared_tolerance = 1E-09;
  sqpnp::PnPSolver solver(*_3dpoints_, *_projections_, wghts, params);

  if(solver.IsValid()){
    solver.Solve();

    poses->resize(solver.NumberOfSolutions());
    for (int i = 0; i < solver.NumberOfSolutions(); i++)
    {
      const sqpnp::SQPSolution *sol = solver.SolutionPtr(i);

//      std::cout << "\nSolution " << i << ":\n";
//      std::cout << *sol << std::endl;
//      std::cout << " Average squared projection error : " << solver.AverageSquaredProjectionErrors().at(i) << std::endl;

      (*poses)[i].block<3,3>(0,0) = Eigen::Map<const Eigen::Matrix<double, 3, 3, Eigen::RowMajor> >( sol->r_hat.data() );
      (*poses)[i].block<3,1>(0,3) = Eigen::Map<const Eigen::Matrix<double, 3, 1, Eigen::ColMajor> >( sol->t.data() );
    }
  }
  else return 0;

  return solver.NumberOfSolutions();
}
#endif

// minimal solver with SQPnP
int PoseEstimator::MinimalSolver(const std::vector<int>& sample,
                                 std::vector<robust_pose_pnp::Matrix34d, Eigen::aligned_allocator<robust_pose_pnp::Matrix34d> > *poses) const {

  //const int npts = sample.size();
  //std::cout<< "in minimal solver (sqpnp) "<<npts<<"\n"<<std::flush;

  poses->clear();

#ifdef ROBUST_POSE_PNP_DEBUG
  if (npts < samplesz_) return 0;
#endif

  return PoseEstimator::sqpnp_solve(sample, poses);
}


// non minimal solver with SQPnP
int PoseEstimator::NonMinimalSolver(const std::vector<int>& sample,
                                    robust_pose_pnp::Matrix34d *pose) const {
  const int npts = sample.size();
//std::cout<< "in non minimal solver "<<npts<<"\n"<<std::flush;

  // note that RansacLib might call the non minimal solver with fewer than nm_samplesz_ points, hence the following check
  if (npts < nm_samplesz_) return 0;

  std::vector<Matrix34d, Eigen::aligned_allocator<Matrix34d> > poses;

  int n = PoseEstimator::sqpnp_solve(sample, &poses);

  if (n > 1) std::cerr << "More than one solution in PoseEstimator::NonMinimalSolver()!\n" << std::flush;

  *pose = poses[0];

  return n > 0;
}

// Evaluates the pose on the i-th point pair.
double PoseEstimator::EvaluateModelOnPoint(const robust_pose_pnp::Matrix34d& pose,
                                           int i) const {
  const sqpnp::_Point& xyz = (*_3dpoints_)[i];

#if 0
  double Xc     =       pose(0, 0)*xyz.vector[0] + pose(0, 1)*xyz.vector[1] + pose(0, 2)*xyz.vector[2] + pose(0, 3),
         Yc     =       pose(1, 0)*xyz.vector[0] + pose(1, 1)*xyz.vector[1] + pose(1, 2)*xyz.vector[2] + pose(1, 3),
         inv_Zc = 1./ ( pose(2, 0)*xyz.vector[0] + pose(2, 1)*xyz.vector[1] + pose(2, 2)*xyz.vector[2] + pose(2, 3) );

#else
  Eigen::Vector3d prod = pose.block<3,3>(0,0) * xyz.vector + pose.block<3,1>(0,3);
  double Xc     =     prod(0),
         Yc     =     prod(1),
         inv_Zc = 1./ prod(2);


  double dx = Xc*inv_Zc - (*_projections_)[i].vector[0];
  double dy = Yc*inv_Zc - (*_projections_)[i].vector[1];
#endif

  return dx*dx + dy*dy;
}

}  // namespace robust_pose_pnp
