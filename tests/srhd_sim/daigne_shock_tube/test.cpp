/*
  Test run
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "spatial_distribution.hpp"
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "spatial_reconstruction.hpp"
#include "rigid_wall.hpp"
#include "pcm.hpp"
#include "van_leer.hpp"
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
    eos_(5./3.),
    rs_(eos_.getAdiabaticIndex()),
    bc_(rs_),
    sr_(),
    geometry_(),
    sim_(linspace(0,1,200),
	 Uniform(1),
	 Step(1000,0.1,0.5),
	 Uniform(0),
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
    const RigidWall<simple_vector> bc_;
    PCM<simple_vector, simple_vector> sr_;
    const Planar geometry_;
    SRHDSimulation<simple_vector, simple_vector> sim_;
  };
}

int main()
{
#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL
  SimData sim_data;
  auto& sim = sim_data.getSim();

  main_loop
    (sim,
     SafeTimeTermination<simple_vector,simple_vector>(0.303,1e6),
     &SRHDSimulation<simple_vector,simple_vector>::timeAdvance,
     WriteTime<simple_vector,simple_vector>("time.txt"));

#ifdef PARALLEL
  write_hdf5_snapshot(sim,"final_"+int2str(get_mpi_rank())+".h5");
#else
  write_hdf5_snapshot<simple_vector, simple_vector>(sim,"final.h5");
#endif // PARALLEL

  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
