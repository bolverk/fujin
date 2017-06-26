/*
  Relativistic shock tube
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <fenv.h>
#include "srhd_sim.hpp"
#include "van_leer.hpp"
#include "ideal_gas.hpp"
#include "hll.hpp"
#include "rigid_wall.hpp"

namespace {
void WritePrimitives(SRHDSimulation const& sim, string const& fname)
{
  std::ofstream f(fname.c_str());
  for(size_t i=0;i<sim.getHydroSnapshot().cells.size();++i){
    const Primitive p = sim.getHydroSnapshot().cells[i];
    f << sim.GetCellCentre(i) << " ";
    f << p.Density << " ";
    f << p.Pressure << " ";
    f << celerity2velocity(p.Celerity) << std::endl;
  }
  f.close();
}

void WriteMidVals(SRHDSimulation const& sim, string const& fname)
{
  std::ofstream f(fname.c_str());
  const Primitive p = sim.getHydroSnapshot().cells[50];
  f << p.Pressure << "\n";
  f << celerity2velocity(p.Celerity) << "\n";
  f.close();
}
}

using namespace std;

int main()
{
  feenableexcept(FE_INVALID | 
		 FE_DIVBYZERO | 
		 FE_OVERFLOW  | 
		 FE_UNDERFLOW);

  vector<double> vertex;
  const size_t n = 100;
  vertex.resize(n);
  for (size_t i=0; i<n; i++)
    vertex[i] = (1.0*i)/(1.0*n-1.0);

  double g = 4./3.;
  Uniform dd(1.0);
  Step dp(0.2,0.1,0.5); 
  Uniform dv(0.0); 
  IdealGas eos(g);
  HLL rs(eos);
  VanLeer sr;
  RigidWall bc(rs);
  const Planar geometry;

  SRHDSimulation sim(vertex,
		     dd, dp, dv,
		     bc, bc,
		     eos,
		     rs,
		     sr,
		     geometry);

  // Main process
  size_t iter = 50;
  vector<double> e(iter);
  while(sim.GetTime()<0.7){
    sim.TimeAdvance();
  }

  // Write data to file
  WriteMidVals(sim,"res.txt");
  WritePrimitives(sim,"plot.txt");

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
