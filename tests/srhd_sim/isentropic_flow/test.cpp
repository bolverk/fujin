/*
  Isentropic flow - verifies the Riemann invariant are properly advected
  See following article for more detail
  W. Zhang & A. I. MacFadyen
  "RAM: A Relativistic Adaptive mesh refinement Hydrodynamics Code"
  Astrophys. Journ. Supp. Ser. 164:255-279 (2006)
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "spatial_distribution.hpp"
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "pcm.hpp"
#include "rigid_wall.hpp"
#include "collela.hpp"
#include "isentropic_flow.hpp"
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
  vector<double> vertex = linspace(-0.35,1,200);
  double g = 5./3.;
  Collela dd2(1,1,0.3,0);
  ConstEntropy dp2(100,g,dd2);
  ConstRiemannInv dv2
    (calc_riemann_invariant(g,dd2(-10),dp2(-10),0,-1),
     g, dd2, dp2);
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  PCM<simple_vector, simple_vector> sr;
  RigidWall bc(rs);
  const double tf = 0.8;
  Planar planar;

  SRHDSimulation<simple_vector, simple_vector> sim
    (vertex,
     dd2, dp2, dv2,
     bc, bc,
     eos,
     rs,
     sr,
     planar);

#ifdef PARALLEL
  write_hdf5_snapshot(sim, "initial_"+int2str(get_mpi_rank())+".h5");
#else
  write_hdf5_snapshot(sim, "initial.h5");
#endif // PARALLEL

  // Main process
  while(sim.getTime()<tf){
    sim.timeAdvance();
  }

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
