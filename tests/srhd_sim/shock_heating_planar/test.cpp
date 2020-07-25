/*
  Planar shock heating
  See following article for more detail
  W. Zhang, A. MacFadyen, "RAM: A Relativistic Adaptive mesh refinement Hydrodynamics Code", 
  The Astrophysical Journal Supplement Series, 164:255-279, 2006 May.
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "spatial_distribution.hpp"
#include "ideal_gas.hpp"
#include "linear_rs.hpp"
#include "pcm.hpp"
#include "van_leer.hpp"
#include "rigid_wall.hpp"
#include "free_flow.hpp"
#include "diagnostics.hpp"
#include <fenv.h>
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

using namespace std;

class SimData
{
public:

  SimData(void):
    dd_(1), 
    dp_(1e-6), 
    dv_(velocity2celerity(-(1.-1e-9))),
    eos_(5./3.),
    rs_(eos_), 
    left_bc_(rs_),
    right_bc_(),
    sr_(),
    geometry_(),
    sim_(linspace(0,1,100),
	 dd_,dp_,dv_,
	 left_bc_, right_bc_,
	 eos_,rs_,sr_,
	 geometry_) {}

  SRHDSimulation& GetSim(void)
  {
    return sim_;
  }

private:

  Uniform dd_;
  Uniform dp_;
  Uniform dv_;
  IdealGas eos_;
  LinearRS rs_;
  RigidWall left_bc_;
  FreeFlow right_bc_;
  //PCM sr_;
  VanLeer sr_;
  const Planar geometry_;
  SRHDSimulation sim_;
};

namespace {

  void main_loop(SRHDSimulation& sim)
  {
    double tf = 0.5;

    while(sim.GetTime()<tf){

      sim.TimeAdvance();
    }
  }

  void Output(SRHDSimulation const& sim)
  {
    ofstream f("res.txt");
    f << sim.getHydroSnapshot().cells[sim.getHydroSnapshot().cells.size()/2].Pressure << "\n";
    f.close();
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

  SimData simd;
  SRHDSimulation& sim = simd.GetSim();

  main_loop(sim);

#ifdef PARALLEL
  if(get_mpi_rank()==0)
#endif // PARALLEL
  Output(sim);

  // Finalise
#ifdef PARALLEL
  if(get_mpi_rank()==0)
#endif // PARALLEL
  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
