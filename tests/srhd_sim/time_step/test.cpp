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
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

//#define N 100
//template<class T> using CE = array<T, N+1>;
//template<class T> using CP = array<T, N>;
//template<class T> using CE = vector<T>;
//template<class T> using CP = vector<T>;

namespace {
  void WritecalcTimeStep
  (const SRHDSimulation<simple_vector, simple_vector>
   & sim,
   string const& fname)
  {
    std::ofstream f(fname.c_str());
    f << sim.calcTimeStep() << "\n";
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
  Step dp(0.2,0.1,0.5);
  Uniform dv(0.0);
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  PCM<simple_vector, simple_vector> sr;
  RigidWall<simple_vector> bc(rs);
  const Spherical geometry;

  SRHDSimulation<simple_vector, simple_vector>sim
    (vertex,
     dd, dp, dv,
     bc, bc,
     eos,
     rs,
     sr,
     geometry);

  // Write data to file
#ifdef PARALLEL
  WritecalcTimeStep(sim, "res_"+int2str(get_mpi_rank())+".txt");
#else
  WritecalcTimeStep(sim,"res.txt");
#endif // PARALLEL

  // Finalise
  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
