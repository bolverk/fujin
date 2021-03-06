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

#define N 1000
template<class T> using CE = array<T, N+1>;
template<class T> using CP = array<T, N>;

using namespace std;

class SimData
{
public:

  SimData(void):
    //    dd_(1), 
    //    dp_(1e-6), 
    //    dv_(velocity2celerity(-(1.-1e-9))),
    eos_(5./3.),
    rs_(eos_), 
    left_bc_(rs_),
    right_bc_(),
    sr_(),
    geometry_(),
    sim_(linspace<N+1>(0,1),
	 [](double){return 1;},
	 [](double){return 1e-6;},
	 [](double){return velocity2celerity(-(1-1e-9));},
	 //	 dd_,dp_,dv_,
	 left_bc_, right_bc_,
	 eos_,rs_,sr_,
	 geometry_) {}

  auto& GetSim(void)
  {
    return sim_;
  }

private:

  //  Uniform dd_;
  //  Uniform dp_;
  //  Uniform dv_;
  IdealGas eos_;
  LinearRS rs_;
  RigidWall<CP> left_bc_;
  FreeFlow<CP> right_bc_;
  //PCM sr_;
  VanLeer<CE, CP> sr_;
  const Planar geometry_;
  SRHDSimulation<CE, CP> sim_;
};

namespace {

  void main_loop(SRHDSimulation<CE, CP>& sim)
  {
    double tf = 0.5;

    while(sim.getTime()<tf){

      sim.timeAdvance();
    }
  }

  void Output
  (const SRHDSimulation<CE, CP>& sim)
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
  auto& sim = simd.GetSim();

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
