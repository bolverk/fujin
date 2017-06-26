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

using namespace std;

namespace {
void WritePrimitives(SRHDSimulation const& sim, string const& fname)
{
  ofstream f(fname.c_str());
  string dlmtr = " "; // Delimiter
  for(size_t i=0;i<sim.getHydroSnapshot().cells.size();++i){
    f << sim.GetCellCentre(i) << dlmtr;
    Primitive p = sim.getHydroSnapshot().cells[i];
    f << p.Density << dlmtr;
    f << p.Pressure << dlmtr;
    f << celerity2velocity(p.Celerity) << endl;
  }
  f.close();
}

void WriteMidVals(SRHDSimulation const& sim, string const& fname)
{
  ofstream f;
  f.open(fname.c_str());
  Primitive p = sim.getHydroSnapshot().cells[sim.getHydroSnapshot().cells.size()/2+1];
  f << p.Pressure << endl;
  f << celerity2velocity(p.Celerity) << endl;
  f.close();
}
}

int main()
{
  const int n = 200;
  vector<double> vertex = linspace(0,1,n);

  double g = 5./3.;
  Uniform dd(1.0);
  Step dp(1000.0,0.1,0.5);
  Uniform dv(0.0);
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  VanLeer sr;
  const Planar geometry;
  RigidWall bc(rs);

  SRHDSimulation sim(vertex,
		     dd, dp, dv,
		     bc, bc,
		     eos,
		     rs,
		     sr,
		     geometry);

  // Main process
  double tf= 0.4;
  while(sim.GetTime()<tf){
    sim.TimeAdvance2ndOrder();    
  }

  // Write data to file
  WriteMidVals(sim,"res.txt");
  WritePrimitives(sim, "plot.txt");

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
