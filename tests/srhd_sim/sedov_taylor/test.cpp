/*
  Sedov Taylor explosion
*/

#include <cassert>
#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "spatial_distribution.hpp"
#include "ideal_gas.hpp"
//#include "imgrs.hpp"
#include "hll.hpp"
#include "pcm.hpp"
#include "rigid_wall.hpp"
#include "main_loop.hpp"
#include <fenv.h>

using namespace std;

namespace {
  double get_position_max_pressure(SRHDSimulation const& sim)
  {
    double pmax = sim.getHydroSnapshot().cells[0].Pressure;
    double rpmax = sim.GetCellCentre(0);
    for(size_t i=0;i<sim.getHydroSnapshot().cells.size();i++){
      Primitive p = sim.getHydroSnapshot().cells[i];
      if(pmax<p.Pressure){
	pmax = p.Pressure;
	rpmax = sim.GetCellCentre(i);
      }
    }
    return rpmax;
  }

  void WriteVectors2File(vector<double> const& v1,
			 vector<double> const& v2,
			 string const& fname)
  {
    ofstream f;
    f.open(fname.c_str());
    for(size_t i=0;i<v1.size();i++){
      f << v1[i] << " " << v2[i] << endl;
    }
    f.close();
  }

  class SimData
  {
  public:

    SimData(void):
      eos_(5./3.),
      rs_(eos_),
      bc_(rs_),
      sr_(),
      geometry_(),
      sim_(linspace(0,1,1000),
	   Uniform(1),
	   Step(1e-3, 1e-9, 0.02),
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
    IdealGas eos_;
    HLL rs_;
    RigidWall bc_;
    PCM sr_;
    const Spherical geometry_;
    SRHDSimulation sim_;
  };

  class ShockFrontTracker: public DiagnosticFunction
  {
  public:
    
    ShockFrontTracker(const string& fname):
      fname_(fname), positions_(), times_() {}

    void operator()(const SRHDSimulation& sim) const
    {
      times_.push_back(sim.GetTime());
      positions_.push_back(get_position_max_pressure(sim));
    }

    ~ShockFrontTracker(void)
    {
      WriteVectors2File(times_, positions_, fname_);
    }
    
  private:
    const string fname_;
    mutable vector<double> positions_;
    mutable vector<double> times_;
  };

  void my_main_loop(SRHDSimulation& sim)
  {
    ShockFrontTracker sft("rpmax.txt");
    SafeTimeTermination stt(50, 1e6);
    main_loop(sim, stt, &SRHDSimulation::TimeAdvance, sft);
  }
}

int main()
{
  feenableexcept(FE_INVALID   | 
		 FE_DIVBYZERO | 
		 FE_OVERFLOW  | 
		 FE_UNDERFLOW);

  my_main_loop(SimData().getSim());

  ofstream("test_terminated_normally.res").close();
  return 0;
}
