/*
  Checks that the program calculates cross section areas properly
  in spherical geometry
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "spatial_distribution.hpp"
#include "ideal_gas.hpp"
#include "imgrs.hpp"
#include "utilities.hpp"
#include "pcm.hpp"
#include "rigid_wall.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

using namespace std;

int main()
{

#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL

  vector<double> vertex;
  const size_t n = 1000;
  vertex.resize(n);
  for (size_t i=0; i<n; i++)
    vertex[i] = 1.0+static_cast<double>(i)/static_cast<double>(n);

  double g = 4./3.;
  Uniform dd(1.0);
  Step dp(0.1,0.2,0.9);
  Uniform dv(0.0);
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  PCM<simple_vector, simple_vector> sr;
  RigidWall<simple_vector> bc(rs);
  const Spherical geometry;

  SRHDSimulation<simple_vector, simple_vector> sim
    (vertex,
     dd, dp, dv,
     bc,
     bc,
     eos,
     rs,
     sr,
     geometry);

  // Write data to file
#ifdef PARALLEL
  if(get_mpi_rank()==0){
#endif // PARALLEL
    ofstream f;
    f.open("res.txt");
    double vol1 = 0;
    double vol2 = 0;
    for (size_t i=0;i<sim.getHydroSnapshot().cells.size();i++)
      {
	vol1 = sim.getVolume(i+1)-sim.getVolume(i);
	vol2 = sim.getArea(i)*
	  (sim.getHydroSnapshot().edges[i+1]-
	   sim.getHydroSnapshot().edges[i]);
	f<<vol1<<" "<<vol2<<endl;
      }
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
