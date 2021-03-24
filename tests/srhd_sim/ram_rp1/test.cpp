/*
  RAM test number 1.
  See following article for more detail
  W. Zhang, A. MacFadyen, "RAM: A Relativistic Adaptive mesh refinement Hydrodynamics Code", 
  The Astrophysical Journal Supplement Series, 164:255-279, 2006 May.
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "pcm.hpp"
#include "ideal_gas.hpp"
#include "hll.hpp"
#include "rigid_wall.hpp"
#include "hdf5_snapshot.hpp"
#include "universal_error.hpp"
#include "main_loop.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

#if SCAFFOLDING != 1
using CE = vector<double>;
using CP = vector<Primitive>;
#endif // SCAFFOLDING

using namespace std;

namespace {

  void main_loop(SRHDSimulation
#if SCAFFOLDING != 1
		 <CE, CP>
#endif // SCAFFOLDING
		 & sim)
  {
    SafeTimeTermination
#if SCAFFOLDING != 1
      <CE, CP>
#endif // SCAFFOLDING
      term_cond(0.4, 1e5);
    WriteTime
      #if SCAFFOLDING != 1
      <CE, CP>
#endif // SCAFFOLDING
      diag("time.txt");
    main_loop(sim,
	      term_cond,
	      &SRHDSimulation
#if SCAFFOLDING != 1
	      <CE, CP>
#endif // SCAFFOLDING
	      ::timeAdvance,
	      diag);
#ifdef PARALLEL
    write_hdf5_snapshot(sim,"final_"+int2str(get_mpi_rank())+".h5");
#else
    write_hdf5_snapshot(sim,"final.h5");
#endif // PARALLEL
  }
}

class SimData
{
public:

  SimData(void):
    adiabatic_index_(5./3.),
    edges_(linspace(0,1,100)),
    density_(10,1,0.5),
    pressure_(13.33,1e-8,0.5),
    velocity_(0),
    eos_(adiabatic_index_),
    rs_(eos_),
    sr_(),
    bc_(rs_),
    geometry_(),
    sim_(edges_,
	 density_,
	 pressure_,
	 velocity_,
	 bc_, bc_,
	 eos_, rs_,
	 sr_, 
	 geometry_) {}

  auto& getSim(void)
  {
    return sim_;
  }

private:
  double adiabatic_index_;
  vector<double> edges_;
  Step density_;
  Step pressure_;
  Uniform velocity_;
  IdealGas eos_;
  //  IdealGasRiemannSolver rs_;
  HLL rs_;
  PCM sr_;
  RigidWall bc_;
  const Planar geometry_;
  SRHDSimulation
#if SCAFFOLDING != 1
  <CE, CP>
#endif // SCAFFOLDING
  sim_;
};

namespace {

  void main_loop_weh(SRHDSimulation
#if SCAFFOLDING != 1
		     <CE, CP>
#endif // SCAFFOLDING
		     & sim)
  {
    try{
      main_loop(sim);
    }
    catch(UniversalError const& eo){
      report_error(eo,cout);
      throw;
    }
  }

}

int main()
{
#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL
  SimData sim_data;
  auto& sim = sim_data.getSim();

  main_loop_weh(sim);

  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
