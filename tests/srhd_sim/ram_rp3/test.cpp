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
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

using namespace std;

int main()
{
#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL

  vector<double> vertex;
  const size_t n = 100;
  vertex.resize(n);
  for (size_t i=0; i<n; i++)
    vertex[i] = static_cast<double>(i)/static_cast<double>(n-1);

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

  SRHDSimulation<simple_vector, simple_vector> sim
    (vertex,
     dd, dp, dv,
     LeftBC,
     RightBC,
     eos,
     rs,
     sr,
     geometry);

  // Main process
  while(sim.getTime()<tf){
    sim.timeAdvance();
  }

  // Output
#ifdef PARALLEL
  write_hdf5_snapshot(sim,"final_"+int2str(get_mpi_rank())+".h5");
#else
  write_hdf5_snapshot(sim,"final.h5");
#endif // PARALLEL

  // Finalise
  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
