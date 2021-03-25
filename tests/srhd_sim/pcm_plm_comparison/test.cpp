/*
  Relativistic shock tube
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "spatial_distribution.hpp"
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

using CE = vector<double>;
using CP = vector<Primitive>;

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
  Step dp(0.2,0.1,0.5); 
  Uniform dv(0.0); 
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  HLL rs2(eos);
  PCM sr;
  VanLeer sr2;
  RigidWall bc(rs);
  double tf = 0.8;
  const Planar geometry;

  SRHDSimulation<CE, CP> sim
    (vertex,
     dd, dp, dv,
     bc, bc,
     eos,
     rs,
     sr,
     geometry);
  SRHDSimulation<CE, CP> sim2
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
  write_hdf5_snapshot<CE, CP>(sim, "pcm.h5");
  write_hdf5_snapshot<CE, CP>(sim, "plm.h5");

  // Finalise
  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
