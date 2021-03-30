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
#include "parallel_helper.hpp"
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
	 //	 Uniform(1),
	 [](double){return 1;},
	 //	 Step(0.2,0.1,0.5),
	 [](double x){return x>0.5?0.1:0.2;},
	 //	 Uniform(0),
	 [](double){return 0;},
	 bc_,
	 bc_,
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
    const Planar geometry_;
    const RigidWall<simple_vector> bc_;
    SRHDSimulation<simple_vector,simple_vector> sim_;
  };
}

using namespace std;

int main()
{
#ifdef PARALLEL
  MPI_Init(NULL,NULL);
#endif // PARALLEL

  SimData sim_data;
  auto& sim = sim_data.getSim();

  // Main process
  main_loop(sim,
	    IterationTermination<simple_vector, simple_vector>(10),
	    &SRHDSimulation<simple_vector,simple_vector>::timeAdvance,
#ifdef PARALLEL
	    TotalEnergyHistory("res_"+int2str(get_mpi_rank())+".txt"));
#else

  TotalEnergyHistory<simple_vector, simple_vector>("res.txt"));

#endif // PARALLEL

  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

 return 0;
}
