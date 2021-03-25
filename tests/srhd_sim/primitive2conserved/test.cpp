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
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

using CE = vector<double>;
using CP = vector<Primitive>;

using namespace std;

int main()
{

#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL

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

  SRHDSimulation<CE, CP> sim
    (vertex,
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
#ifdef PARALLEL
  if(get_mpi_rank()==0){
#endif // PARALLEL
    ofstream f;
    f.open("res.txt");
    f << cv.Mass << endl;
    f << cv.Momentum << endl;
    f << cv.Energy << endl;
    f.close();
#ifdef PARALLEL
  }
#endif // PARALLEL

  // Finalise
  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
