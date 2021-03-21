#include <iostream>
#include <cmath>
#include "spatial_reconstruction.hpp"
#include "spatial_distribution.hpp"
#include "periodic.hpp"
#include "ideal_gas.hpp"
#include "hll.hpp"
#include "van_leer.hpp"
#include "utilities.hpp"
#include "srhd_sim.hpp"
#include "diagnostics.hpp"
#include "main_loop.hpp"
#include "hdf5_snapshot.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

using namespace std;

namespace {

  class SimData
  {
  public:

    SimData(void):
      eos_(4./3.),
      rs_(eos_),
      bc_(rs_),
      sr_(),
      geometry_(),
      sim_(linspace(0,1,100),
	   SineWave(1,1,0,2),
	   Uniform(1),
	   Uniform(0.5),
	   bc_, bc_,
	   eos_,
	   rs_,
	   sr_,
	   geometry_) {}

    SRHDSimulation& getSim(void)
    {
      return sim_;
    }

  private:
    const IdealGas eos_;
    const HLL rs_;
    const Periodic bc_;
    VanLeer sr_;
    const Planar geometry_;
    SRHDSimulation sim_;
  };
}

int main(void)
{
#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL
  SimData sim_data;
  SRHDSimulation& sim = sim_data.getSim();

  //  write_snapshot(sim,"init_cond.txt");
#ifdef PARALLEL
  write_hdf5_snapshot(sim, "initial_"+int2str(get_mpi_rank())+".h5");
#else
  write_hdf5_snapshot(sim, "initial.h5");
#endif // PARALLEL

  main_loop(sim,
	    SafeTimeTermination(2,1e6),
	    &SRHDSimulation::timeAdvance,
	    WriteTime("time.txt"));

  //  write_snapshot(sim,"snapshot.txt");
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
