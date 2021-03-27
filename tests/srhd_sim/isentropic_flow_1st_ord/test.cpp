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
#include "hdf5_snapshot.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

using namespace std;

namespace {
  class density_distribution: public SpatialDistribution
  {
  private:
    double dref;
    double alpha;
    double l;
  public:
    density_distribution(double idref, double ialpha, double il):
      dref(idref),alpha(ialpha),l(il) {}
    double operator()(double x) const
    {
      double f = 0;
      if(abs(x)<l)
	f = pow(pow(x/l,2)-1,4);
      return dref*(1+alpha*f);
    }
  };

  class pressure_distribution: public SpatialDistribution
  {
  private:
    double k;
    double g;
    const density_distribution& rdd;
  public:
    pressure_distribution
    (double ik, double ig, const density_distribution& irdd):
      k(ik),g(ig), rdd(irdd) {}
    double operator()(double x) const
    {
      double d = rdd(x);
      return k*pow(d,g);
    }
  };

  double energy(double g, double d, double p)
  {
    return d+p/(g-1);
  }

  double soundspeed(double g, double d, double p)
  {

    double e = energy(g, d, p);
    return sqrt(g*p/(e+p));
  }

  double riemann_invariant(double g, double d, double p,
			   double v, double s)
  {
    double cs = soundspeed(g, d, p);
    return atanh(v)+s*(2/sqrt(g-1))*
      atanh(cs/sqrt(g-1));
  }


  class velocity_distribution: public SpatialDistribution
  {
  private:
    double jm;
    double g;
    const density_distribution& rdd;
    const pressure_distribution& rpd;
  public:
    velocity_distribution
    (double ijm, double ig, const density_distribution& irdd, 
     const pressure_distribution& irpd):
      jm(ijm), g(ig), rdd(irdd), rpd(irpd) {}
  
    double operator()(double x) const
    {
      double d = rdd(x);
      double p = rpd(x);
      double aux = jm - riemann_invariant(g, d, p, 0, -1);
      return tanh(aux);
    }
  };
}

int main()
{
#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL
  vector<double> vertex;
  const size_t n = 200;
  vertex.resize(n);
  for (size_t i=0; i<n; i++)
    vertex[i] = -0.35 + 1.35*static_cast<double>(i)/static_cast<double>(n-1);

  double g = 5./3.;
  density_distribution dd(1,1,0.3);
  pressure_distribution dp(100, g, dd);
  velocity_distribution dv
    (riemann_invariant(g,dd(-10),dp(-10),0,-1),
     g, dd, dp);
  IdealGas eos(g);
  IdealGasRiemannSolver rs(g);
  PCM<simple_vector, simple_vector> sr;
  RigidWall<simple_vector> bc(rs);
  const double tf = 0.8;
  Planar geometry;

  SRHDSimulation<simple_vector, simple_vector> sim
    (vertex,
     dd, dp, dv,
     bc, bc,
     eos,
     rs,
     sr,
     geometry);

#ifdef PARALLEL
  write_hdf5_snapshot(sim, "initial_"+int2str(get_mpi_rank())+".h5");
#else
  write_hdf5_snapshot(sim, "initial.h5");
#endif // PARALLEL

  // Main process
  while(sim.getTime()<tf){
    sim.timeAdvance1stOrder();
  }

#ifdef PARALLEL
  write_hdf5_snapshot(sim, "final_"+int2str(get_mpi_rank())+".h5");
#else
  write_hdf5_snapshot(sim, "final.h5");
#endif // PARALLEL

  // Finalise
  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
