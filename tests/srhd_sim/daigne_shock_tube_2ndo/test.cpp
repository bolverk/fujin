#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "spatial_distribution.hpp"
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "spatial_reconstruction.hpp"
#include "rigid_wall.hpp"
#include "pcm.hpp"
#include "van_leer.hpp"
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

  const int n = 200;
  vector<double> vertex = linspace(0,1,n);

  double g = 5./3.;
  Uniform dd(1.0);
  Step dp(1000.0,0.1,0.5);
  Uniform dv(0.0);
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  VanLeer sr;
  const Planar geometry;
  RigidWall bc(rs);

  SRHDSimulation sim(vertex,
		     dd, dp, dv,
		     bc, bc,
		     eos,
		     rs,
		     sr,
		     geometry);

  // Main process
  double tf = 0.4;
  while(sim.GetTime()<tf){
    sim.TimeAdvance();    
  }

  // Write data to file
#ifdef PARALLEL
  write_hdf5_snapshot(sim, "final_"+int2str(get_mpi_rank())+".h5");
#else
  write_hdf5_snapshot(sim, "final.h5");
#endif // PARALLEL

  // Finalise
  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
