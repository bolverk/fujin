/*
  Checks that the file does compile
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "pcm.hpp"
#include "rigid_wall.hpp"
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
  //  Uniform dd(1.0);
  const auto dd = [](double){return 1;};
  //  Step dp(0.2,0.1,0.5);
  const auto dp = [](double x){return x<0.5?0.2:0.1;};
  //  Uniform dv(0.0);
  const auto dv = [](double){return 0;};
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  PCM<simple_vector, simple_vector> sr;
  RigidWall<simple_vector> bc(rs);
  Spherical geometry;

  SRHDSimulation<simple_vector, simple_vector> sim
    (vertex,
     dd, dp, dv,
     bc, bc,
     eos,
     rs,
     sr,
     geometry);
  sim.timeAdvance();

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
