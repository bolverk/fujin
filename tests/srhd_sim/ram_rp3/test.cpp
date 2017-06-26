/*
  RAM test number 3.
  See following article for more detail
  W. Zhang, A. MacFadyen, "RAM: A Relativistic Adaptive mesh refinement Hydrodynamics Code", 
  The Astrophysical Journal Supplement Series, 164:255-279, 2006 May.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "van_leer.hpp"
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "spatial_reconstruction.hpp"
#include "rigid_wall.hpp"
#include "free_flow.hpp"
#include "hdf5_snapshot.hpp"

using namespace std;

int main()
{
  vector<double> vertex;
  const size_t n = 100;
  vertex.resize(n);
  for (size_t i=0; i<n; i++)
    vertex[i] = (1.0*i)/(1.0*n-1.0);

  double g = 4./3.;
  Uniform dd(1.0);
  Step dp(1.0,10.0,0.5);
  Step dv(velocity2celerity(0.9),0.0,0.5);
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  VanLeer sr;
  FreeFlow LeftBC;
  RigidWall RightBC(rs);
  const double tf = 0.3;
  Planar geometry;

  SRHDSimulation sim(vertex,
		     dd, dp, dv,
		     LeftBC,
		     RightBC,
		     eos,
		     rs,
		     sr,
		     geometry);

  // Main process
  while(sim.GetTime()<tf){
    sim.TimeAdvance();
  }

  // Output
  write_hdf5_snapshot(sim,"final.h5");

  // Finalise
  ofstream("test_terminated_normally.res").close();
  return 0;
}
