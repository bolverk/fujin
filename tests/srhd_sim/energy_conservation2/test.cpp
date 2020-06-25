/*
  Checks that the file does compile
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "pcm.hpp"
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "spatial_reconstruction.hpp"
#include "rigid_wall.hpp"
#include "diagnostics.hpp"

namespace {
void WriteVector(vector<double> const& v, string const& fname)
{
  std::ofstream f(fname.c_str());
  for(size_t i=0;i<v.size();i++){
    f << v[i] << "\n";
  }
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
    vertex[i] = static_cast<double>(i)/static_cast<double>(n-1);

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
		     bc,
		     bc,
		     eos,
		     rs,
		     sr,
		     geometry);

  // Main process
  size_t iter = 10;
  vector<double> e(iter);
  for (size_t i=0;i<iter;i++)
    {
      sim.TimeAdvance();
      e[i] = TotalEnergy(sim);
    }

  // Write data to file
  WriteVector(e,"res.txt");
  WritePrimitives(sim,"plot.txt");

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
