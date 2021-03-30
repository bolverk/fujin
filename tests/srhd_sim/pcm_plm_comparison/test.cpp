/*
  Relativistic shock tube
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "hll.hpp"
#include "pcm.hpp"
#include "van_leer.hpp"
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
  auto dd = [](double){return 1;};
  auto dp = [](double x){return x<0.5?0.2:0.1;};
  auto dv = [](double){return 0;};
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  HLL rs2(eos);
  PCM<simple_vector, simple_vector> sr;
  VanLeer<simple_vector, simple_vector> sr2;
  RigidWall<simple_vector> bc(rs);
  double tf = 0.8;
  const Planar geometry;

  SRHDSimulation<simple_vector, simple_vector> sim
    (vertex,
     dd, dp, dv,
     bc, bc,
     eos,
     rs,
     sr,
     geometry);
  SRHDSimulation<simple_vector, simple_vector> sim2
    (vertex,
      dd, dp, dv,
      bc, bc,
      eos,
      rs2,
      sr2,
      geometry);

  // Main process
  while(sim.getTime()<tf)
    sim.timeAdvance();

  while(sim2.getTime()<tf)
    sim2.timeAdvance();

  // Write data to file
  write_hdf5_snapshot<simple_vector, simple_vector>(sim, "pcm.h5");
  write_hdf5_snapshot<simple_vector, simple_vector>(sim, "plm.h5");

  // Finalise
  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
