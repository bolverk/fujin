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
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "pcm.hpp"
#include "rigid_wall.hpp"
#include "isentropic_flow.hpp"
#include "hdf5_snapshot.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

using namespace std;

namespace {
  double calc_riemann_invariant
  (double g, double d, double p,
   double v, double s)
  {
    const double cs = IdealGas(g).dp2ba(d,p);
    return atanh(v)+s*(2/sqrt(g-1))*
      atanh(cs/sqrt(g-1));
  }
}

int main()
{
#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL
  vector<double> vertex = linspace(-0.35,1,200);
  double g = 5./3.;
  auto dd2 = [](double x){
	       const double f = abs(x)<0.3?pow(pow(x/0.3,2)-1,4):0;
	       return 1+f;
	     };
  auto dp2 = [&dd2,&g](double x){return 100*pow(dd2(x),g);};
  const double jm = calc_riemann_invariant
    (g, dd2(-10), dp2(-10), 0, -1);
  auto dv2 = [&g, &dd2, &dp2, &jm](double x){
	       const double d = dd2(x);
	       const double p = dp2(x);
	       const double aux = jm - calc_riemann_invariant(g,d,p,0,-1);
	       return atanh(aux);
	     };
  //  Collela dd2(1,1,0.3,0);
  //  ConstEntropy dp2(100,g,dd2);
  /*
    ConstRiemannInv dv2
    (calc_riemann_invariant(g,dd2(-10),dp2(-10),0,-1),
    g, dd2, dp2);
  */
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  PCM<simple_vector, simple_vector> sr;
  RigidWall<simple_vector> bc(rs);
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
