/*
  RAM test number 2.
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
#include "spatial_reconstruction.hpp"
#include "rigid_wall.hpp"
#include "hdf5_snapshot.hpp"
#include "main_loop.hpp"

using namespace std;

namespace {
  class SimData
  {
  public:

    SimData(void):
      eos_(5./3.),
      rs_(eos_),
      bc_(rs_),
      geometry_(),
      sim_(linspace(0,1,100),
	   Uniform(1),
	   Step(1e3,1e-2,0.5),
	   Uniform(0),
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
    const RigidWall bc_;
    PCM sr_;
    const Planar geometry_;
    SRHDSimulation sim_;
  };

  void my_main_loop(SRHDSimulation& sim)
  {
    SafeTimeTermination term_cond(0.4, 1e6);
    WriteTime diag("time.txt");
    main_loop(sim,
	      term_cond,
	      &SRHDSimulation::TimeAdvance,
	      diag);
    write_hdf5_snapshot(sim,"final.h5");
  }
}

int main()
{
  my_main_loop(SimData().getSim());

  ofstream("test_terminated_normally.res");
  return 0;
}
