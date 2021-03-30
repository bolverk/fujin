#include <iostream>
#include <cmath>
#include "spatial_reconstruction.hpp"
#include "periodic.hpp"
#include "ideal_gas.hpp"
#include "linear_rs.hpp"
#include "van_leer.hpp"
#include "utilities.hpp"
#include "srhd_sim.hpp"
#include "diagnostics.hpp"
#include "pcm.hpp"
#include "main_loop.hpp"
#include "hdf5_snapshot.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

template<class T> using simple_vector = vector<T>;

using namespace std;

namespace {

  double calc_riemann_invariant
  (double g, double d, double p,
   double v, double s)
  {
    const double cs = IdealGas(g).dp2ba(d,p);
    return atanh(v)+s*(2/sqrt(g-1))*
      atanh(cs/sqrt(g-1));
  }

    class InitCond
  {
  public:

    InitCond(double mean_density,
	     double pert_density,
	     double mean_pressure,
	     double g):
      cd_(mean_density),
      amp_(pert_density),
      cp_(mean_pressure),
      g_(g),
      density_([this](double x){return cd_+amp_*sin(2*M_PI*x);}),
      pressure_([this](double x){return cp_*pow(density_(x)/cd_,g_);}),
      jm_(calc_riemann_invariant
	  (g,
	   density_(0),
	   pressure_(0),
	   0,
	   -1)),
      velocity_([this](double x){
		  const double d = density_(x);
		  const double p = pressure_(x);
		  const double aux = jm_-calc_riemann_invariant(g_,d,p,0,-1);
		  return tanh(aux);}){}

    function<double(double)> const& getDist(string const& dname) const
    {
      if("density"==dname)
	return density_;
      else if("pressure"==dname)
	return pressure_;
      else if("velocity"==dname)
	return velocity_;
      else
	throw "Unknown distribution name";
    }

  private:
    const double cd_;
    const double amp_;
    const double cp_;
    const double g_;
    const function<double(double)> density_;
    const function<double(double)> pressure_;
    const double jm_;
    const function<double(double)> velocity_;
  };

  class SimData
  {
  public:

    SimData(void):
    eos_(4./3.),
    init_cond_
    (1,1e-5,1e-3,
     eos_.getAdiabaticIndex()),
    rs_(eos_),
    bc_(rs_),
    sr_(),
    geometry_(),
    sim_(linspace(0,1,20),
	 init_cond_.getDist("density"),
	 init_cond_.getDist("pressure"),
	 init_cond_.getDist("velocity"),
	 bc_, bc_,
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
    const InitCond init_cond_;
    const LinearRS rs_;
    const Periodic<simple_vector> bc_;
    VanLeer<simple_vector, simple_vector> sr_;
    const Planar geometry_;
    SRHDSimulation<simple_vector, simple_vector> sim_;
  };
}

int main(void)
{
#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL
  SimData sim_data;
  auto& sim = sim_data.getSim();

#ifdef PARALLEL
  write_hdf5_snapshot(sim, "initial_"+int2str(get_mpi_rank())+".h5");
#else
  write_hdf5_snapshot(sim, "initial.h5");
#endif // PARALLEL

  main_loop
    (sim,
     SafeTimeTermination<simple_vector, simple_vector>(10,1e6),
     &SRHDSimulation<simple_vector, simple_vector>::timeAdvance,
     WriteTime<simple_vector, simple_vector>("time.txt"));

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
