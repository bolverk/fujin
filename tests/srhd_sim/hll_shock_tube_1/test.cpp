/*
  Relativistic shock tube
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "spatial_distribution.hpp"
#include "ideal_gas.hpp"
#include "hll.hpp"
#include "pcm.hpp"
#include "rigid_wall.hpp"
#include "main_loop.hpp"
#include "hdf5_snapshot.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

using CE = vector<double>;
using CP = vector<Primitive>;

namespace {
  /*
  void WriteMidVals(SRHDSimulation const& sim, string fname)
  {
    std::ofstream f(fname.c_str());
    const Primitive p = sim.getHydroSnapshot().cells[sim.getHydroSnapshot().cells.size()/2];
    f << p.Pressure << "\n";
    f << celerity2velocity(p.Celerity) << "\n";
    f.close();
  }
  */

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
	   Uniform(1),
	   Step(0.2,0.1,0.5),
	   Uniform(0),
	   bc_, bc_,
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
    const HLL rs_;
    const RigidWall bc_;
    PCM sr_;
    const Planar geometry_;
    SRHDSimulation<CE, CP> sim_;
  };
}

using namespace std;

int main()
{
#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL
  SimData sim_data;
  auto& sim = sim_data.getSim();

    main_loop(sim,
	      IterationTermination<CE,CP>(50),
	      &SRHDSimulation<CE, CP>::timeAdvance,
	      WriteTime<CE,CP>("time.txt"));

  // Output
  //WriteMidVals(sim,"res.txt");
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
