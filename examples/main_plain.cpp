//
// main.cpp
//
// Manolis Lourakis (lourakis **at** ics forth gr), September 2020
// 
// Example demo program for SQPnP with data points loaded from plain 2D arrays
//

#include <vector>
#include <iostream>
#include <chrono>

#include <types.h>
#include <sqpnp.h>


int main()
{
  const int n=7;

  double pts3[n][3]={
  {-0.429857595273321, -0.441798127281825, 0.714342354521372},
  {-2.1568268264648, 0.113521604867983, -0.148634122716948},
  {0.694636908485644, -0.737067927134015, -1.38877746946909},
  {-1.07051455287146, -1.2122304801284, -0.841002964233812},
  {0.509844073252947, -1.07097319594739, 0.675410167109412},
  {0.40951585099, 2.2300713816052, 0.365229861025625},
  {2.04320214188098, 1.11847674401846, 0.623432173763436},
  };

  double pts2[n][2]={
#if 0
  // no noise
  {0.139024436737141, -0.00108631784422283},
  {0.149897105048989, 0.270584578309815},
  {-0.118448642309468, -0.0844116551810971},
  {0.0917181969674735, 0.0435196877212059},
  {0.100243308685939, -0.178506520365217},
  {-0.296312157121094, 0.220675975198136},
  {-0.331509880499455, -0.213091587841007},
#else
  // noisy
  {0.138854772853285, -0.00157437083896972},
  {0.149353089173631, 0.269826809256435},
  {-0.118391028248405, -0.0842834292914752},
  {0.0937833539430025, 0.0473371294380393},
  {0.101410594775151, -0.179030803711188},
  {-0.294749181228375, 0.221134043355639},
  {-0.334084299358372, -0.21071853326318},
#endif
  };

  auto start = std::chrono::high_resolution_clock::now();

  std::vector<sqpnp::_Point> _3dpoints(n);
  std::vector<sqpnp::_Projection> _projections(n);

  for (int i = 0; i < n; i++) {
    const double *p3=pts3[i];
    const double *p2=pts2[i];

    _3dpoints[i]=sqpnp::_Point(p3[0], p3[1], p3[2]);
    _projections[i]=sqpnp::_Projection(p2[0], p2[1]);
  }

  // demonstration of passing parameters to the solver
  sqpnp::SolverParameters params;
  params.omega_nullspace_method = sqpnp::OmegaNullspaceMethod::RRQR;
  // equal weights for all points
  sqpnp::PnPSolver solver(_3dpoints, _projections, std::vector<double>(n, 1.0), params);

  auto stop = std::chrono::high_resolution_clock::now();

  if (solver.IsValid()) {
    solver.Solve();
    stop = std::chrono::high_resolution_clock::now();
    std::cout << "SQPnP found " << solver.NumberOfSolutions() << " solution(s)"<< std::endl;
    for (int i = 0; i < solver.NumberOfSolutions(); i++)
    {
      std::cout << "\nSolution " << i << ":\n";
      std::cout << *solver.SolutionPtr(i) << std::endl;
      std::cout << " Average squared projection error : " << solver.AverageSquaredProjectionErrors()[i] << std::endl;
    }
  }
  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
  std::cout << "Time taken by SQPnP: " << duration.count() << " microseconds" << std::endl << std::endl;

  return 0;
}
