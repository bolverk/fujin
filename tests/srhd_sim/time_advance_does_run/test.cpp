/*
  Checks that the simulation can complete a single cycle
 */

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "pcm.hpp"
#include "rigid_wall.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

#define N 100
template<class T> using CE = array<T, N+1>;
template<class T> using CP = array<T, N>;

using namespace std;

int main()
{

#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL

  CE<double> vertex;
  const size_t n = N+1;
  //  vertex.resize(n);
  for (size_t i=0; i<n; i++)
    vertex[i] = static_cast<double>(i)/static_cast<double>(n-1);

  double g = 4./3.;
  //  Uniform dd(1.0);
  const auto dd = [](double){return 1;};
  //  Step dp(0.1,0.2,0.5);
  const auto dp = [](double x){return x<0.5?0.1:0.2;};
  //  Uniform dv(0.0);
  const auto dv = [](double){return 0;};  
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  PCM<CE, CP> sr;
  RigidWall<CP> bc(rs);
  const Spherical geometry;

  SRHDSimulation<CE, CP>
    sim(vertex,
		     dd, dp, dv,
		     bc, bc,
		     eos,
		     rs,
		     sr,
		     geometry);
  sim.timeAdvance();

  // Write data to file
#ifdef PARALLEL
  if(get_mpi_rank()==0){
#endif // PARALLEL
  ofstream f;
  f.open("res.txt");
  f << 0;
  f.close();
#ifdef PARALLEL
  }
#endif // PARALLEL

  // Finalise
  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

 return 0;
}
