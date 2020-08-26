#include <iostream>
#include <fstream>
#include <vector>
#include "srhd_sim.hpp"
#include "pcm.hpp"
#include "ideal_gas.hpp"
#include "hll.hpp"
#include "rigid_wall.hpp"
#include "hdf5_snapshot.hpp"
#include "universal_error.hpp"
#include "main_loop.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

using namespace std;

namespace {

  void main_loop(SRHDSimulation& sim,
		 const int& inum,
		 const string& outname)
  {
    //    SafeTimeTermination term_cond(tf, 1e5);
    IterationTermination term_cond(inum);
    WriteTime diag("time.txt");
    main_loop(sim,
	      term_cond,
	      &SRHDSimulation::TimeAdvance,
	      diag);
#ifdef PARALLEL
    write_hdf5_snapshot(sim,"final_"+int2str(get_mpi_rank())+".h5");
#else
    write_hdf5_snapshot(sim,outname);
#endif // PARALLEL
  }
}

class SimData
{
public:

  SimData(void):
    adiabatic_index_(5./3.),
    edges_(linspace(0,1,100)),
    density_(10,1,0.5),
    pressure_(13.33,1e-8,0.5),
    velocity_(0),
    eos_(adiabatic_index_),
    rs_(eos_),
    sr_(),
    bc_(rs_),
    geometry_(),
    sim_(edges_,
	 density_,
	 pressure_,
	 velocity_,
	 bc_, bc_,
	 eos_, rs_,
	 sr_, 
	 geometry_) {}

  SimData(const string& fname):
    adiabatic_index_(5./3.),
    edges_(linspace(0,1,100)),
    density_(10,1,0.5),
    pressure_(13.33,1e-8,0.5),
    velocity_(0),
    eos_(adiabatic_index_),
    rs_(eos_),
    sr_(),
    bc_(rs_),
    geometry_(),
    sim_(read_hdf5_snapshot(fname),
	 bc_, bc_,
	 eos_, rs_,
	 sr_, 
	 geometry_) {}

  SRHDSimulation& getSim(void)
  {
    return sim_;
  }

private:
  double adiabatic_index_;
  vector<double> edges_;
  Step density_;
  Step pressure_;
  Uniform velocity_;
  IdealGas eos_;
  //  IdealGasRiemannSolver rs_;
  HLL rs_;
  PCM sr_;
  RigidWall bc_;
  const Planar geometry_;
  SRHDSimulation sim_;
};

namespace {

  void main_loop_weh(SRHDSimulation& sim,
		     const int& inum,
		     const string& outname)
  {
    try{
      main_loop(sim, inum, outname);
    }
    catch(UniversalError const& eo){
      report_error(eo,cout);
      throw;
    }
  }

}

int main()
{
#ifdef PARALLEL
  MPI_Init(NULL, NULL);
  throw("not implemented yet");
#endif // PARALLEL
  {
    SimData sim_data;
    SRHDSimulation& sim = sim_data.getSim();

    main_loop_weh(sim, 300, "continuous.h5");
  }

  {
    SimData sim_data;
    SRHDSimulation& sim = sim_data.getSim();

    main_loop_weh(sim, 150, "checkpoint.h5");
  }

  {
    SimData sim_data("checkpoint.h5");
    SRHDSimulation& sim = sim_data.getSim();

    main_loop_weh(sim, 149, "interrupted.h5");
  }
  
  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL

  return 0;
}
