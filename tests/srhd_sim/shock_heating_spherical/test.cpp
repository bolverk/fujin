/*
  Spherical shock heating
  See following article for more detail
  F. Daigne, R. Mochkovitch
  "Gamma-ray bursts from internal shocks in a relativistic wind:
  a hydrodynamic study"
  Atron. Astrophys. 358, 1157-1166 (2000)
 */

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "spatial_distribution.hpp"
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "pcm.hpp"
#include "rigid_wall.hpp"
#include "free_flow.hpp"

namespace {
void WritePrimitives(SRHDSimulation const& sim, string const& fname)
{
  std::ofstream f(fname.c_str());
  string dlmtr = " "; // Delimiter
  for(size_t i=0;i<sim.getHydroSnapshot().cells.size();i++){
    f << sim.GetCellCentre(i) << dlmtr;
    const Primitive p = sim.getHydroSnapshot().cells[i];
    f << p.Density << dlmtr;
    f << p.Pressure << dlmtr;
    f << celerity2velocity(p.Celerity) << std::endl;
  }
  f.close();
}

  /*
void WriteConserveds(SRHDSimulation const& sim, string const& fname)
{
  std::ofstream f(fname.c_str());
  for(size_t i=0;i<sim.getHydroSnapshot().cells.size();i++){
    Conserved c = sim.GetConserved(i);
    f << c.Mass << " ";
    f << c.Momentum << " ";
    f << c.Energy << "\n";
  }
  f.close();
}
  */
}

using namespace std;

int main()
{
  vector<double> vertex;
  const size_t n = 200;
  vertex.resize(n);
  for (size_t i=0; i<n; i++)
    vertex[i] = (1.0*i)/(1.0*n-1.0);

  double g = 4./3.;
  const double vin = 0.99999;
  const double win = velocity2celerity(vin);
  const double rsr = 0.01;
  Uniform dd(1.0);
  Uniform dp(1e-6);
  Step dv(0,-win,rsr); 
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  PCM sr;
  RigidWall LeftBC(rs);
  FreeFlow RightBC;
  const double tf = 0.5;
  Spherical geometry;

  SRHDSimulation sim(vertex,
		     dd, dp, dv,
		     LeftBC, RightBC,
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
  f << sim.getHydroSnapshot().cells[10].Pressure << "\n";
  f.close();

  WritePrimitives(sim, "plot.txt");
  //  WriteConserveds(sim, "plot_rs.txt");

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
