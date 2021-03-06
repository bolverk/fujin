/*
  Checks that the program calculates the volumes properly
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
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
  const size_t n = 100;
  vertex.resize(n);
  for (size_t i=0; i<n; i++)
    vertex[i] = static_cast<double>(i)/static_cast<double>(n);

  double g = 4./3.;
  /*
  Uniform dd(1.0);
  Step dp(0.1,0.2,0.5);
  Uniform dv(0.0);
  */
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  PCM<simple_vector, simple_vector> sr;
  RigidWall<simple_vector> bc(rs);
  const Spherical geometry;

  SRHDSimulation <simple_vector, simple_vector> sim
    (vertex,
     [](double){return 1;},
     [](double x){return x<0.5?0.1:0.2;},
     [](double){return 0;},
     //     dd, dp, dv,
     bc, bc,
     eos,
     rs,
     sr,
     geometry);

  // Write data to file
  ofstream f;
#ifdef PARALLEL
  f.open("res_"+int2str(get_mpi_rank())+".txt");
#else
  f.open("res.txt");
#endif // PARALLEL
  //  double vol = 0;
  for (size_t i=0;i<sim.getHydroSnapshot().cells.size();i++)
    {
      const double vol = (4.0*M_PI/3.0)*(pow(sim.getHydroSnapshot().edges[i],3));
      f<<sim.getVolume(i)<<" "<<vol<<endl;
    }
  f.close();

  // Finalise
  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
