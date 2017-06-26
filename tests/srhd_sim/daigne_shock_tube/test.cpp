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

using namespace std;

namespace {
  void WriteMidVals(SRHDSimulation const& sim, string const& fname)
  {
    ofstream f(fname.c_str());
    const Primitive p = sim.getHydroSnapshot().cells[sim.getHydroSnapshot().cells.size()/2];
    f << p.Pressure << endl;
    f << celerity2velocity(p.Celerity) << endl;
    f.close();
  }

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

    SRHDSimulation& getSim(void)
    {
      return sim_;
    }

  private:
    const IdealGas eos_;
    const IdealGasRiemannSolver rs_;
    const RigidWall bc_;
    PCM sr_;
    const Planar geometry_;
    SRHDSimulation sim_;
  };
}

int main()
{
  SimData sim_data;
  SRHDSimulation& sim = sim_data.getSim();

  main_loop(sim,
	    SafeTimeTermination(0.303,1e6),
	    &SRHDSimulation::TimeAdvance,
	    WriteTime("time.txt"));

  write_hdf5_snapshot(sim,"final.h5");
  WriteMidVals(sim,"res.txt");

  ofstream("test_terminated_normally.res").close();
 return 0;
}
