/*
  Checks that the file does compile
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "pcm.hpp"
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "spatial_reconstruction.hpp"
#include "rigid_wall.hpp"
#include "diagnostics.hpp"
#include "main_loop.hpp"
#ifdef PARALLEL
#include "mpi.h"
#endif // PARALLEL

namespace {

  class SimData
  {
  public:

    SimData(void):
    eos_(4./3.),
    rs_(eos_.getAdiabaticIndex()),
    sr_(),
    geometry_(),
    bc_(rs_),
    sim_(linspace(0,1,100),
	 Uniform(1),
	 Step(0.2,0.1,0.5),
	 Uniform(0),
	 bc_,
	 bc_,
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
    const IdealGasRiemannSolver rs_;
    PCM sr_;
    const Planar geometry_;
    const RigidWall bc_;
    SRHDSimulation sim_;
  };
}

using namespace std;

int main()
{
#ifdef PARALLEL
  MPI_Init(NULL,NULL);
#endif // PARALLEL

  SimData sim_data;
  SRHDSimulation& sim = sim_data.getSim();

  // Main process
  main_loop(sim,
	    IterationTermination(10),
	    &SRHDSimulation::TimeAdvance,
	    TotalEnergyHistory("res.txt"));

  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

 return 0;
}
