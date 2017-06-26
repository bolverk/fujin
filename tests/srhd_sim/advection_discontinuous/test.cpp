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
	 Step(2,1,0.5),
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
  SimData sim_data;
  SRHDSimulation& sim = sim_data.getSim();
  write_snapshot(sim,"init_cond.txt");

  main_loop(sim,
	    SafeTimeTermination(2,1e6),
	    &SRHDSimulation::TimeAdvance,
	    WriteTime("time.txt"));

  write_snapshot(sim,"snapshot.txt");

  ofstream("test_terminated_normally.res").close();
 return 0;
}
