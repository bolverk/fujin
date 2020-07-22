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
#include "advanced_hydrodynamic_variables.hpp"

using namespace std;

int main()
{

  vector<double> vertex;
  const size_t n = 100;
  vertex.resize(n+1);
  for (size_t i=0; i<n; i++)
    vertex[i] = static_cast<double>(i)/static_cast<double>(n);

  double g = 4./3.;
  Uniform dd(1.0);
  Step dp(0.1,0.2,0.5);
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

  const Primitive hs(1.0,0.1, velocity2celerity(0.2));
  Conserved cv = Primitive2Conserved(hs, eos);

  // Write data to file
  ofstream f;
  f.open("res.txt");
  f << cv.Mass << endl;
  f << cv.Momentum << endl;
  f << cv.Energy << endl;
  f.close();

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
