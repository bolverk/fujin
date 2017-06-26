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

namespace {
  void WriteMidVals(SRHDSimulation const& sim, string const& fname)
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
    rs_(eos_.getAdiabaticIndex()),
    sr_(),
    bc_(rs_),
    geometry_(),
    sim_(linspace(0,1,100),
	 Uniform(1),
	 Step(0.2,0.1,0.5),
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
    PCM sr_;
    const RigidWall bc_;
    const Planar geometry_;
    SRHDSimulation sim_;
  };
}

using namespace std;

int main()
{
  feenableexcept(FE_INVALID   | 
		 FE_DIVBYZERO | 
		 FE_OVERFLOW  | 
		 FE_UNDERFLOW);

  SimData sim_data;
  SRHDSimulation& sim = sim_data.getSim();

  main_loop(sim,
	    SafeTimeTermination(0.7,1e6),
	    &SRHDSimulation::TimeAdvance,
	    WriteTime("time.txt"));

  WriteMidVals(sim,"res.txt");
  write_snapshot(sim,"final.txt",14);

  ofstream("test_terminated_normally.res").close();
 return 0;
}
