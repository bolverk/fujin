/*
  Checks that the file does compile
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

namespace {
void WriteTimeStep(SRHDSimulation const& sim, string const& fname)
{
  std::ofstream f(fname.c_str());
  f << sim.TimeStep() << "\n";
  f.close();
}

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
}

using namespace std;

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
  PCM sr;
  RigidWall bc(rs);
  const Spherical geometry;

  SRHDSimulation sim(vertex,
		     dd, dp, dv,
		     bc, bc,
		     eos,
		     rs,
		     sr,
		     geometry);

  // Write data to file
  WriteTimeStep(sim,"res.txt");
  WritePrimitives(sim,"plot.txt");

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
