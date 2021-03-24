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

#if SCAFFOLDING != 1
using CE = vector<double>;
using CP = vector<Primitive>;
#endif // SCAFFOLDING

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
    const RigidWall bc_;
    PCM sr_;
    const Planar geometry_;
    SRHDSimulation
#if SCAFFOLDING != 1
    <CE, CP>
#endif // SCAFFOLDING
    sim_;
  };
}

int main()
{
#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL
  SimData sim_data;
  auto& sim = sim_data.getSim();

#if SCAFFOLDING == 1
  main_loop(sim,
	    SafeTimeTermination(0.303,1e6),
	    &SRHDSimulation::timeAdvance,
	    WriteTime("time.txt"));
#else
  main_loop
    (sim,
     SafeTimeTermination<CE,CP>(0.303,1e6),
     &SRHDSimulation<CE,CP>::timeAdvance,
     WriteTime<CE,CP>("time.txt"));
#endif // SCAFFOLDING

#ifdef PARALLEL
  write_hdf5_snapshot(sim,"final_"+int2str(get_mpi_rank())+".h5");
#else
  write_hdf5_snapshot
#if SCAFFOLDING != 1
    <CE, CP>
#endif // SCAFFOLDING
    (sim,"final.h5");
#endif // PARALLEL

  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
