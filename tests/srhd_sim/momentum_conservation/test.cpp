/*
  Checks that the file does compile
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
#include "diagnostics.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

#if SCAFFOLDING != 1
using CE = vector<double>;
using CP = vector<Primitive>;
#endif // SCAFFOLDING

namespace {
  void WriteVector(vector<double> const& v, string const& fname)
  {
    std::ofstream f(fname.c_str());
    f.precision(14);
    for(size_t i=0;i<v.size();i++){
      f << v[i] << "\n";
    }
    f.close();
  }
}

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
  TwoSteps dp(0.1,0.4,0.2,0.6,0.1);
  TwoSteps dv(0.0,0.4,0.1,0.6,0.0);
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  PCM sr;
  RigidWall bc(rs);
  const Planar geometry;

  SRHDSimulation
#if SCAFFOLDING != 1
    <CE, CP>
#endif // SCAFFOLDING
    sim(vertex,
		     dd, dp, dv,
		     bc, bc,
		     eos,
		     rs,
		     sr,
		     geometry);

  // Main process
  size_t iter = 10;
  vector<double> e(iter);
  for (size_t i=0;i<iter;i++)
    {
      sim.timeAdvance();
      e[i] = TotalMomentum(sim);
    }

  // Write data to file
#ifdef PARALLEL
  WriteVector(e,"res_"+int2str(get_mpi_rank())+".txt");
#else
  WriteVector(e,"res.txt");
#endif // PARALLEL

  // Finalise
  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
