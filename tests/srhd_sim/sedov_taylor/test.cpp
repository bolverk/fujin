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
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

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

  double get_max_pressure(const SRHDSimulation& sim)
  {
    double pmax = sim.getHydroSnapshot().cells[0].Pressure;
    for(size_t i=0;i<sim.getHydroSnapshot().cells.size();++i)
      pmax = max(pmax, sim.getHydroSnapshot().cells[i].Pressure);
    return pmax;
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
      pressures_.push_back(get_max_pressure(sim));
      positions_.push_back(get_position_max_pressure(sim));
    }

    ~ShockFrontTracker(void)
    {
      ofstream f(fname_.c_str());
      for(size_t i=0;i<positions_.size();++i)
	f << times_.at(i) << " "
	  << positions_.at(i) << " "
	  << pressures_.at(i) << endl;
      f.close();
    }
    
  private:
    const string fname_;
    mutable vector<double> positions_;
    mutable vector<double> pressures_;
    mutable vector<double> times_;
  };

  void my_main_loop(SRHDSimulation& sim)
  {
#ifdef PARALLEL
    ShockFrontTracker sft("rpmax_"+int2str(get_mpi_rank())+".txt");
#else
    ShockFrontTracker sft("rpmax.txt");
#endif // PARALLEL
    SafeTimeTermination stt(50, 1e6);
    main_loop(sim, stt, &SRHDSimulation::TimeAdvance, sft);
  }
}

int main()
{
#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL
  feenableexcept(FE_INVALID   | 
		 FE_DIVBYZERO | 
		 FE_OVERFLOW  | 
		 FE_UNDERFLOW);

  my_main_loop(SimData().getSim());

  ofstream("test_terminated_normally.res").close();
  
#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
