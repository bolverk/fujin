/*
  Isentropic flow - verifies the Riemann invariant are properly advected
  See following article for more detail
  W. Zhang & A. I. MacFadyen
  "RAM: A Relativistic Adaptive mesh refinement Hydrodynamics Code"
  Astrophys. Journ. Supp. Ser. 164:255-279 (2006)
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
#include "collela.hpp"
#include "isentropic_flow.hpp"

using namespace std;

namespace {
void write_snapshot(SRHDSimulation const& sim, string const& fname)
{
  ofstream f;
  f.open(fname.c_str());
  for(size_t i=0;i<sim.getHydroSnapshot().cells.size();++i){
    Primitive p = sim.getHydroSnapshot().cells[i];
    f << sim.GetCellCentre(i) << " ";
    f << p.Density << " ";
    f << p.Pressure << " ";
    f << celerity2velocity(p.Celerity) << endl;
  }
  f.close();
}
}

int main()
{
  vector<double> vertex = linspace(-0.35,1,200);
  double g = 5./3.;
  Collela dd2(1,1,0.3,0);
  ConstEntropy dp2(100,g,dd2);
  ConstRiemannInv dv2
    (calc_riemann_invariant(g,dd2(-10),dp2(-10),0,-1),
     g, dd2, dp2);
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  PCM sr;
  RigidWall bc(rs);
  const double tf = 0.8;
  Planar planar;

  SRHDSimulation sim(vertex,
		     dd2, dp2, dv2,
		     bc, bc,
		     eos,
		     rs,
		     sr,
		     planar);

  write_snapshot(sim, "plot0.txt");

  // Main process
  while(sim.GetTime()<tf){
    sim.TimeAdvance();
  }

  write_snapshot(sim, "plot.txt");

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
