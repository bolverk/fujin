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

namespace {
  void WriteMidVals(SRHDSimulation const& sim, string fname)
  {
    std::ofstream f(fname.c_str());
    const Primitive p = sim.getHydroSnapshot().cells[sim.getHydroSnapshot().cells.size()/2];
    f << p.Pressure << "\n";
    f << celerity2velocity(p.Celerity) << "\n";
    f.close();
  }

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
}

using namespace std;

int main()
{
  SimData sim_data;
  SRHDSimulation& sim = sim_data.getSim();

  main_loop(sim,
	    IterationTermination(50),
	    &SRHDSimulation::TimeAdvance,
	    WriteTime("time.txt"));

  // Output
  WriteMidVals(sim,"res.txt");

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
