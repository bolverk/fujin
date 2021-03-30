/*
  Relativistic shock tube
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "pcm.hpp"
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "rigid_wall.hpp"
#include "main_loop.hpp"
#include "diagnostics.hpp"
#include <fenv.h>
#include "hdf5_snapshot.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL_HELPER

namespace {

  class SimData
  {
  public:

    SimData(void):
      eos_(4./3.),
      rs_(eos_.getAdiabaticIndex()),
      sr_(),
      bc_(rs_),
      geometry_(),
      sim_(linspace(0,1,100),
	   //	   Uniform(1),
	   [](double){return 1;},
	   //	   Step(0.2,0.1,0.5),
	   [](double x){return x<0.5?0.2:0.1;},
	   //	   Uniform(0),
	   [](double){return 0;},
	   bc_,bc_,
	   eos_,
	   rs_,
	   sr_,
	   geometry_) {}

    auto& getSim(void)
    {
      return sim_;
    }

  private:
    const IdealGas eos_;
    const IdealGasRiemannSolver rs_;
    PCM<simple_vector, simple_vector> sr_;
    const RigidWall<simple_vector> bc_;
    const Planar geometry_;
    SRHDSimulation<simple_vector, simple_vector> sim_;
  };
}

using namespace std;

int main()
{

#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL

  feenableexcept(FE_INVALID   | 
		 FE_DIVBYZERO | 
		 FE_OVERFLOW  | 
		 FE_UNDERFLOW);

  SimData sim_data;
  auto& sim = sim_data.getSim();

  main_loop(sim,
	    SafeTimeTermination<simple_vector, simple_vector>(0.7,1e6),
	    &SRHDSimulation<simple_vector, simple_vector>::timeAdvance,
	    WriteTime<simple_vector, simple_vector>("time.txt"));

#ifdef PARALLEL
  write_hdf5_snapshot(sim, "final_"+int2str(get_mpi_rank())+".h5");
#else
  write_hdf5_snapshot(sim, "final.h5");
#endif // PARALLEL

  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
