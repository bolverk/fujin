/*
  Relativistic shock tube
 */

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "spatial_distribution.hpp"
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "hll.hpp"
#include "pcm.hpp"
#include "van_leer.hpp"
#include "rigid_wall.hpp"

using namespace std;

namespace {
void WritePrimitives(SRHDSimulation const& sim, string const& fname)
{
  ofstream f;
  f.open(fname.c_str());
  for(size_t i=0;i<sim.getHydroSnapshot().cells.size();++i){
    const Primitive p = sim.getHydroSnapshot().cells[i];
    f << sim.GetCellCentre(i) << " ";
    f << p.Density << " ";
    f << p.Pressure << " ";
    f << velocity2celerity(p.Celerity) << endl;
  }
  f.close();
}

  /*
void WriteConserveds(SRHDSimulation const& sim, string const& fname)
{
  ofstream f;
  f.open(fname.c_str());
  for(size_t i=0;i<sim.getHydroSnapshot().cells.size();++i){
    Conserved c = sim.GetConserved(i);
    f << c.Mass << " ";
    f << c.Momentum << " ";
    f << c.Energy << endl;
  }
  f.close();
}
  */

void WriteRes(SRHDSimulation const& sim1, SRHDSimulation const& sim2)
{
  ofstream f;
  f.open("res.txt");
  f << sim1.getHydroSnapshot().cells[50].Pressure << endl;
  f << velocity2celerity(sim1.getHydroSnapshot().cells[50].Celerity) << endl;
  f << sim2.getHydroSnapshot().cells[50].Pressure << endl;
  f << velocity2celerity(sim2.getHydroSnapshot().cells[50].Celerity) << endl;
  f.close();
}
}

int main()
{
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
  IdealGasRiemannSolver rs(g);
  HLL rs2(eos);
  PCM sr;
  VanLeer sr2;
  RigidWall bc(rs);
  double tf = 0.8;
  const Planar geometry;

  SRHDSimulation sim(vertex,
		     dd, dp, dv,
		     bc, bc,
		     eos,
		     rs,
		     sr,
		     geometry);
  SRHDSimulation sim2(vertex,
		      dd, dp, dv,
		      bc, bc,
		      eos,
		      rs2,
		      sr2,
		      geometry);

  // Main process
  while(sim.GetTime()<tf)
    sim.TimeAdvance();

  while(sim2.GetTime()<tf)
    sim2.TimeAdvance();

  // Write data to file
  WriteRes(sim, sim2);

  WritePrimitives(sim,"plot.txt");
  WritePrimitives(sim2,"plot2.txt");
  //  WriteConserveds(sim,"plot3.txt");
  //  WriteConserveds(sim2,"plot4.txt");

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
