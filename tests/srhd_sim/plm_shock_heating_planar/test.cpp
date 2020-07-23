/*
  Planar shock heating
  See following article for more detail
  W. Zhang, A. MacFadyen, "RAM: A Relativistic Adaptive mesh refinement Hydrodynamics Code", 
  The Astrophysical Journal Supplement Series, 164:255-279, 2006 May.
 */

#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "spatial_distribution.hpp"
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "pcm.hpp"
#include "rigid_wall.hpp"
#include "free_flow.hpp"
#include <fenv.h>
#include "hdf5_snapshot.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

using namespace std;

namespace {
void WritePrimitives(SRHDSimulation const& sim, string const& fname)
{
  ofstream f;
  f.open(fname.c_str());
  for(size_t i=0;i<sim.getHydroSnapshot().cells.size();i++){
    Primitive p = sim.getHydroSnapshot().cells[i];
    f << sim.GetCellCentre(i) << " ";
    f << p.Density << " ";
    f << p.Pressure << " ";
    f << celerity2velocity(p.Celerity) << endl;
  }
  f.close();
}

  /*
void WriteConserveds(SRHDSimulation const& sim, string const& fname)
{
  ofstream f;
  f.open(fname.c_str());
  for(size_t i=0;i<sim.getHydroSnapshot().cells.size();i++){
    Conserved c = sim.GetConserved(i);
    f << c.Mass << " ";
    f << c.Momentum << " ";
    f << c.Energy << endl;
  }
  f.close();
}
  */
}

int main()
{
#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL

  feenableexcept(FE_INVALID | 
		 FE_DIVBYZERO | 
		 FE_OVERFLOW  | 
		 FE_UNDERFLOW);

  /*
  vector<double> vertex;
  const size_t n = 101;
  vertex.resize(n);
  for (size_t i=0; i<n; i++)
    vertex[i] = static_cast<double>(i)/static_cast<double>(n-1);
  */
  const vector<double> vertex = linspace(0,1,101);

  double g = 5./3.;
  //  const double vin = 1.0-1e-10;
  const double vin = 1.0 - 1e-10;
  const double win = velocity2celerity(vin);
  Uniform dd(1.0);
  Uniform dp(1e-6);
  Step dv(win,-win,0.5); 
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  PCM sr;
  FreeFlow LeftBC;
  FreeFlow RightBC;
  const double tf = 0.2;
  const Planar geometry;

  SRHDSimulation sim(vertex,
		     dd, dp, dv,
		     LeftBC,
		     RightBC,
		     eos,
		     rs,
		     sr,
		     geometry);

  // Main process
  while(sim.GetTime()<tf){
    sim.TimeAdvance();
  }

  // Write data to file
  ofstream f;
  f.open("res.txt");
  f << sim.getHydroSnapshot().cells[vertex.size()/2].Pressure << endl;
  f.close();
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
