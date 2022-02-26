//
// robust_pose_pnp.h
//
// Manolis Lourakis (lourakis **at** ics forth gr), February 2022

#ifndef ROBUST_POSE_PNP_H_
#define ROBUST_POSE_PNP_H_

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <random>
#include <vector>

#include <Eigen/Core>

#include <types.h>
#include <sqpnp.h>

namespace robust_pose_pnp {

// pose as [R t]
typedef Eigen::Matrix<double, 3, 4> Matrix34d;

// Interfaces to SQPnP for estimating a 6D pose from 3D-2D matches
class PoseEstimator {
 public:
  PoseEstimator(const std::vector<sqpnp::_Point> *_3dpoints, const std::vector<sqpnp::_Projection> *_projections, int samplesz, int nm_samplesz = 30);

  // Estimates pose with RANSAC and supplied inlier percentage and outlier threshold.
  // If requested, fills in inlier and outlier indices
  int ransacfit(double inlPcent, double outlThresh, robust_pose_pnp::Matrix34d& best_pose, std::vector<int> *idxInliers = nullptr, std::vector<int> *idxOutliers = nullptr) const;

  // sets sample sizes
  inline void set_sample_sizes(int samplesz, int nm_samplesz) {
    samplesz_ = (samplesz < 3)? 3 : samplesz; nm_samplesz_ = (nm_samplesz < 2*samplesz)? 2*samplesz : nm_samplesz;
  }

  inline int min_sample_size() const { return samplesz_; }

  inline int non_minimal_sample_size() const { return nm_samplesz_; }

  inline int num_data() const { return (*_3dpoints_).size(); }

  // Returns the estimated poses in a vector
  int MinimalSolver(const std::vector<int>& sample,
                    std::vector<Matrix34d, Eigen::aligned_allocator<Matrix34d> > *poses) const;

  // Returns 0 if no pose could be estimated and 1 otherwise
  int NonMinimalSolver(const std::vector<int>& sample,
                       Matrix34d *pose) const;

  // Evaluates the error for the i-th data point
  double EvaluateModelOnPoint(const Matrix34d &pose, int i) const;

  // "least squares" solver. Calls NonMinimalSolver
  inline void LeastSquares(const std::vector<int>& sample,
                           Matrix34d *pose) const {
    NonMinimalSolver(sample, pose);
  }

 private:
  // Interface to the SQPnP solver for both the min/non minimal cases
  int sqpnp_solve(const std::vector<int>& sample,
                  std::vector<Matrix34d, Eigen::aligned_allocator<Matrix34d> > *poses) const;

 protected:
  // pointers to vectors holding the 3D & 2D points
  const std::vector<sqpnp::_Point> *_3dpoints_;
  const std::vector<sqpnp::_Projection> *_projections_;

  // sample sizes
  int samplesz_, nm_samplesz_;

};

}  // namespace robust_pose_pnp

#endif  // ROBUST_POSE_PNP_H_
