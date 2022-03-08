//
// main.cpp
//
// Manolis Lourakis (lourakis **at** ics forth gr), February 2022
// 
// Demo program for robust camera pose estimation from 3D - 2D correspondences with SQPnP
// The demo uses SQPnP both as a minimal and as a non-minimal solver in a LO-RANSAC framework.
// 
// The program reads from two files: one with the camera intrinsics K and another containing
// the 3D - 2D corresponding points in separate lines as
// X Y Z x y
// with x y in pixels. To use normalized coordinates for x y, specify the identity matrix as K
//
// For more advanced pose estimation (e.g., a non-linear refinement on the inliers),
// check the posest library at https://users.ics.forth.gr/~lourakis/posest/ 
//

#include <cstdio>
#include <vector>
#include <iostream>
#include <chrono>

#include <types.h>
#include <sqpnp.h>

#include "robust_pose_pnp.h"


#define MAXSTRLEN 1024

/* read matching points from a file and normalize the projections */
static int readMatchingPoints(char *fname, Eigen::Matrix<double, 3, 3>& K1, std::vector<sqpnp::_Projection>& pts2D, std::vector<sqpnp::_Point>& pts3D)
{
register int i;
int ncoords, nmatches;
double X, Y, Z, x, y, nx, ny;
FILE *fp;
char buf[MAXSTRLEN], *ptr;
long int fpos;

  if((fp=fopen(fname, "r"))==NULL){
    fprintf(stderr, "cannot open file %s\n", fname);
    exit(1);
  }

  do{
    fpos=ftell(fp);
    ptr=fgets(buf, MAXSTRLEN-1, fp);
    if(!ptr || ferror(fp)){
      fprintf(stderr, "File %s: error reading line \"%s\"\n", fname, ptr);
      exit(1);
    }
  } while(!feof(fp) && buf[0]=='#'); /* skip comments */

  if(feof(fp)){
    fclose(fp);
    return 0;
  }

  ncoords=sscanf(buf, "%lf%lf%lf%lf%lf", &X, &Y, &Z, &x, &y);
  if(ncoords==5){ /* no lines number */
    for(nmatches=1; !feof(fp); nmatches++){
      i=fscanf(fp, "%*g%*g%*g%*g%*g\n");
      if(ferror(fp)){
        fprintf(stderr, "File %s: error reading point coordinates, line %d\n", fname, nmatches + 1);
        exit(1);
      }
    }

    if(fseek(fp, fpos, SEEK_SET)){ /* rewind right after any comment lines */
      fprintf(stderr, "fseek failed in readMatchingPoints()\n");
      exit(1);
    }
  }
  else{
    sscanf(buf, "%d", &nmatches);
  }

  pts2D.reserve(nmatches);
  pts3D.reserve(nmatches);

  /* read in points and store them */
  for(i=0; !feof(fp); i++){
    ncoords=fscanf(fp, "%lf%lf%lf%lf%lf\n", &X, &Y, &Z, &x, &y);
    if(ncoords==EOF) break;

    if(ncoords!=5){
      fprintf(stderr, "File %s: line %d contains only %d coordinates\n", fname, i + 1, ncoords);
      exit(1);
    }
    if(ferror(fp)){
      fprintf(stderr, "File %s: error reading point coordinates, line %d\n", fname, i + 1);
      exit(1);
    }

    pts3D.emplace_back(X, Y, Z);

    nx=K1(0, 0)*x + K1(0, 1)*y + K1(0, 2);
    ny=K1(1, 0)*x + K1(1, 1)*y + K1(1, 2);
    pts2D.emplace_back(nx, ny);
  }
  fclose(fp);

  if(i!=nmatches){
    fprintf(stderr, "number of actual points in file %s does not agree with that in first line (%d != %d)!\n",
                     fname, i, nmatches);
    exit(1);
  }

  return nmatches;
}

/* reads the 3x3 intrinsic calibration matrix contained in a file */
static void readCalibParams(char *fname, double K[9])
{
  FILE *fp;
  int i, ch=EOF;
  char buf[MAXSTRLEN];

  if((fp=fopen(fname, "r"))==NULL){
    fprintf(stderr, "cannot open file %s, exiting\n", fname);
    exit(1);
  }

  while(!feof(fp) && (ch=fgetc(fp))=='#') /* skip comments */
    fgets(buf, MAXSTRLEN-1, fp);

  if(feof(fp)){
    fclose(fp);
    K[0]=K[1]=K[2]=K[3]=K[4]=K[5]=K[6]=K[7]=K[8]=0.0;
    return;
  }

  ungetc(ch, fp);

  for(i=0; i<3; i++){
    if(fscanf(fp, "%lf%lf%lf\n", K, K+1, K+2)!=3){
      fprintf(stderr, "cannot read three numbers from row %d in file %s, exiting\n", i+1, fname);
      exit(1);
    }
    K+=3;
  }

  fclose(fp);
}


int main(int argc, char *argv[])
{
  std::vector<sqpnp::_Point> pts3D;
  std::vector<sqpnp::_Projection> pts2D;

  int npts;
  char *icalfile, *matchesfile;
  double ical[9]={1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0}; // eye(3)

  std::chrono::time_point<std::chrono::high_resolution_clock> start, stop;

  /* arguments parsing */
  if(argc!=3){
    fprintf(stderr, "Usage: %s <K> <matched points>\n", argv[0]);
    exit(1);
  }

  icalfile=argv[1];
  matchesfile=argv[2];

  if(strncmp(icalfile, "-", 1))
    readCalibParams(icalfile, ical);

  Eigen::Matrix<double, 3, 3> K = Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor> > (ical);
  Eigen::Matrix<double, 3, 3> K1= K.inverse();
  npts=readMatchingPoints(matchesfile, K1, pts2D, pts3D);

  std::cout << argv[0] << ": read " << pts3D.size() << " points"<< std::endl;

#if 0
  for(int i=0; i<npts; ++i){
    printf("%g %g %g  %g %g\n", pts3D[i].vector[0], pts3D[i].vector[1], pts3D[i].vector[2], pts2D[i].vector[0], pts2D[i].vector[1]);
  }
#endif


  // the following will run SQPnP on all input points
#if 0
  start = std::chrono::high_resolution_clock::now();

  sqpnp::SolverParameters params;
  //params.sqp_max_iteration=17;
  //params.rank_tolerance=1E-06;
  params.nearest_rotation_method=sqpnp::NearestRotationMethod::FOAM;
  sqpnp::PnPSolver solver(pts3D, pts2D, std::vector<double>(npts, 1.0), params);

  stop = std::chrono::high_resolution_clock::now();

  if(solver.IsValid()){
    solver.Solve();
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "SQPnP on all points: " << solver.NumberOfSolutions() << " solution(s)"<< std::endl;
    for (int i = 0; i < solver.NumberOfSolutions(); i++)
    {
      std::cout << "\nSolution " << i << ":\n";
      std::cout << *solver.SolutionPtr(i) << std::endl;
      std::cout << " Average squared projection error : " << solver.AverageSquaredProjectionErrors().at(i) << std::endl;
    }
  }
  auto micros = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "Time taken by SQPnP: " << micros.count() << " microseconds" << std::endl << std::endl;
#endif


  start = std::chrono::high_resolution_clock::now();

  robust_pose_pnp::Matrix34d bestRt;
  std::vector<int> outidx;

  robust_pose_pnp::PoseEstimator rpe(&pts3D, &pts2D, 4 /*3*/, 50);
  //rpe.set_sample_sizes(4, 30); // changes sample sizes dynamically
  rpe.ransacfit(0.8, 0.01, bestRt, nullptr, &outidx);

  stop = std::chrono::high_resolution_clock::now();

  std::cout << "\nRANSAC estimate" << std::endl;
  std::cout << bestRt << std::endl;

  std::cout << "\nOutlying pairs" << std::endl;
  for(size_t i=0; i < outidx.size(); i++)
    std::cout << outidx[i] << ' ';
  std::cout << std::endl;

  auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  std::cout << "Time taken by RANSAC: " << millis.count() << " milliseconds" << std::endl << std::endl << std::flush;

  return 0;
}
