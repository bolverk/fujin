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
#include "van_leer.hpp"
#include "rigid_wall.hpp"
#include "universal_error.hpp"
#include "diagnostics.hpp"
#include "main_loop.hpp"
#include "hdf5_snapshot.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

#if SCAFFOLDING != 1
using CE = vector<double>;
using CP = vector<Primitive>;
#endif // SCAFFOLDING

using namespace std;

namespace {
  class density_distribution: public SpatialDistribution
  {
  private:
    const double dref;
    const double alpha;
    const double l;
  public:
    density_distribution(double idref, double ialpha, double il):
      dref(idref),alpha(ialpha),l(il) {}
    double operator()(double x) const
    {
      const double f = pow(pow(x/l,2)-1,4)*(abs(x)<l);
      return dref*(1+alpha*f);
    }
  };

  class pressure_distribution: public SpatialDistribution
  {
  private:
    const double k;
    const double g;
    const density_distribution& rdd;
  public:
    pressure_distribution
    (double ik, double ig, const density_distribution& irdd):
      k(ik),g(ig), rdd(irdd) {}
    double operator()(double x) const
    {
      const double d = rdd(x);
      return k*pow(d,g);
    }
  };

  double riemann_invariant(double g, double d, double p,
			   double v, double s)
  {
    const double cs = IdealGas(g).dp2ba(d,p);
    return atanh(v)+s*(2/sqrt(g-1))*
      atanh(cs/sqrt(g-1));
  }


  class velocity_distribution: public SpatialDistribution
  {
  private:
    const double jm;
    const double g;
    const density_distribution& rdd;
    const pressure_distribution& rpd;
  public:
    velocity_distribution
    (double ijm, double ig, const density_distribution& irdd, 
     const pressure_distribution& irpd):
      jm(ijm), g(ig), rdd(irdd), rpd(irpd) {}
  
    double operator()(double x) const
    {
      const double d = rdd(x);
      const double p = rpd(x);
      const double aux = jm - riemann_invariant(g, d, p, 0, -1);
      return tanh(aux);
    }
  };

  class InitialConditions
  {
  public:

    InitialConditions(double g, double far_away=-10):
      dd_(1,1,0.3),
      dp_(100,g,dd_),
      dv_(riemann_invariant(g,
			    dd_(far_away),
			    dp_(far_away),
			    0,-1),
	  g,dd_,dp_) {}

    density_distribution const& getDensity(void)
    {
      return dd_;
    }

    pressure_distribution const& getPressure(void)
    {
      return dp_;
    }

    velocity_distribution const& getVelocity(void)
    {
      return dv_;
    }

  private:
    const density_distribution dd_;
    const pressure_distribution dp_;
    const velocity_distribution dv_;
  };

  class SimData
  {
  public:

    SimData(void):
    eos_(5./3.),
    rs_(eos_.getAdiabaticIndex()),
    sr_(),
    bc_(rs_),
    geometry_(),
    sim_(linspace(-0.35,1.35,200),
	 InitialConditions(eos_.getAdiabaticIndex()).getDensity(),
	 InitialConditions(eos_.getAdiabaticIndex()).getPressure(),
	 InitialConditions(eos_.getAdiabaticIndex()).getVelocity(),
	 bc_,bc_,
	 eos_,
	 rs_,
	 sr_,
	 geometry_) {}

    auto& getSim(void)
    {
      return sim_;
    }

  private:
    const IdealGas eos_;
    const IdealGasRiemannSolver rs_;
    VanLeer sr_;
    const RigidWall bc_;
    const Planar geometry_;
    SRHDSimulation
#if SCAFFOLDING != 1
    <CE, CP>
#endif // SCAFFOLDING
    sim_;
  };
}

int main()
{
#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif //PARALLEL
  
  SimData sim_data;
  auto& sim = sim_data.getSim();

#ifdef PARALLEL
  write_hdf5_snapshot(sim, "initial_"+int2str(get_mpi_rank())+".h5");
#else
  write_hdf5_snapshot(sim, "initial.h5");
#endif // PARALLEL

  main_loop(sim,
	    SafeTimeTermination
#if SCAFFOLDING != 1
	    <CE, CP>
#endif // SCAFFOLDING
	    (0.8,1e6),
#ifdef PARALLEL
	    &SRHDSimulation::timeAdvance,
#else
	    &SRHDSimulation
	    #if SCAFFOLDING != 1
	    <CE, CP>
#endif // SCAFFOLDING
	    ::timeAdvance2ndOrder,
#endif // PARALLEL
	    WriteTime
#if SCAFFOLDING != 1
	    <CE, CP>
#endif // SCAFFOLDING
	    ("time.txt"));

#ifdef PARALLEL
  write_hdf5_snapshot(sim, "final_"+int2str(get_mpi_rank())+".h5");
#else
  write_hdf5_snapshot(sim, "final.h5");
#endif // PARALLEL

  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

 return 0;
}
