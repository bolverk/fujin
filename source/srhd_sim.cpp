/*! \file srhd_sim.cpp
  \brief Relativistic, Godunov type, lagrangian 1d hydrodynamic simulation
  \author Almog Yalinewich
 */

#include <cmath>
#include <complex>
#include <cassert>
#include "srhd_sim.hpp"
#include "utilities.hpp"
#include "srhydro.hpp"
#include "universal_error.hpp"
#include "advanced_hydrodynamic_variables.hpp"
#ifdef PARALLEL
#include "mpi.h"
#endif // PARALLEL
#include "spdlog/spdlog.h"

using namespace std;

double SRHDSimulation::getArea(size_t i) const
{
  const double r = 0.5*(data_.edges.at(i)+data_.edges.at(i+1));
  return geometry_.calcArea(r);
}

namespace{

#ifdef PARALLEL
  int get_mpi_rank(void)
  {
    int res;
    MPI_Comm_rank(MPI_COMM_WORLD, &res);
    return res;
  }

  int get_mpi_size(void)
  {
    int res;
    MPI_Comm_size(MPI_COMM_WORLD, & res);
    return res;
  }
#endif // PARALLEL

  vector<double> distribute_vertices
    (const vector<double>& vertices)
  {
#ifdef PARALLEL
    spdlog::debug("Inside parallel distribute vertices");
    const size_t cell_num = vertices.size()-1;
    const int rank = get_mpi_rank();
    const int size = get_mpi_size();
    spdlog::debug("cell_num {0}, rank {1}, size {2}",
		  cell_num, rank, size);
    
    vector<int> partition(static_cast<size_t>(size), cell_num/size);
    for(size_t i=0;i<cell_num%size;++i)
      ++partition.at(i);
    
    vector<int> cumpar(partition.size()+1,0);
    for(size_t i=1;i<cumpar.size();++i)
      cumpar.at(i) = cumpar.at(i-1) + partition.at(i-1);

    const size_t low = cumpar.at(rank);
    const size_t high = cumpar.at(rank+1);

    vector<double> res(high-low+1,0);
    for(size_t i=0;i<res.size();++i)
      res.at(i) = vertices.at(low+i);
    spdlog::debug("low {0}, high {1}",
		  res.front(),
		  res.back());
    return res;
#else
    return vertices;
#endif // PARALLEL
  }
}

#ifdef PARALLEL
#endif //PARALLEL

SRHDSimulation::SRHDSimulation 
(const vector<double>& rVertices,
 const SpatialDistribution& DensityDistribution,
 const SpatialDistribution& PressureDistribution,
 const SpatialDistribution& VelocityDistribution,
 const BoundaryCondition& piInnerBC,
 const BoundaryCondition& piOuterBC,
 const EquationOfState& reos,
 const RiemannSolver& rRiemannSolver,
 const SpatialReconstruction& pInterpolationMethod,
 const Geometry& geometry):
  data_(distribute_vertices(rVertices),
	InitCells(distribute_vertices(rVertices),
		  DensityDistribution,
		  PressureDistribution,
		  VelocityDistribution)),
  eos_(reos), 
  psvs_(distribute_vertices(rVertices).size(),
	RiemannSolution()),
  rs_(rRiemannSolver), 
  sr_(pInterpolationMethod),
  cfl_(1./3.), 
  consVars_(primitives_to_new_conserveds(data_.cells,reos)),
  restMass_(serial_generate(RestMassCalculator(data_,geometry))),
  geometry_(geometry),
  time_(0),
  cycle_(0),
  innerBC_(piInnerBC), 
  outerBC_(piOuterBC) {}

SRHDSimulation::SRHDSimulation
(const NewHydroSnapshot<vector<double>, vector<Primitive> >& init_cond,
 const BoundaryCondition& piInnerBC,
 const BoundaryCondition& piOuterBC,
 const EquationOfState& reos,
 const RiemannSolver& rRiemannSolver,
 const SpatialReconstruction& pInterpolationMethod,
 const Geometry& geometry):
  data_(init_cond),
  eos_(reos), 
  psvs_(init_cond.edges.size(),
	RiemannSolution()),
  rs_(rRiemannSolver), 
  sr_(pInterpolationMethod),
  cfl_(1./3.), 
  consVars_(primitives_to_new_conserveds(data_.cells,reos)),
  restMass_(serial_generate(RestMassCalculator(data_,geometry))),
  geometry_(geometry),
  time_(0),
  cycle_(0),
  innerBC_(piInnerBC), 
  outerBC_(piOuterBC) {}

void SRHDSimulation::overrideCFL(double cfl_new)
{
  cfl_ = cfl_new;
}

double SRHDSimulation::calcTimeStep(void) const
{
  return cfl_*MaxTimeStep(data_.edges, data_.cells,eos_);
}

double SRHDSimulation::calcTimeStepForCell(size_t i) const
{
  return MaxTimeStep(data_.edges[i+1]-data_.edges[i],
		     data_.cells[i],eos_);
}

const vector<double>& SRHDSimulation::getRestMasses(void) const
{
  return restMass_;
}

double SRHDSimulation::getVolume(size_t i) const
{
  return geometry_.calcVolume(data_.edges.at(i));
}

void SRHDSimulation::timeAdvance1stOrder(void)
{
#ifdef PARALLEL
  const double dt_candidate = calcTimeStep();
  double temp = dt_candidate;
  MPI_Allreduce(&dt_candidate,&temp,1,
		MPI_DOUBLE,
		MPI_MIN,
		MPI_COMM_WORLD);
  const double dt = temp;
  spdlog::debug("dt = {0}",dt);
#else
  const double dt = calcTimeStep();
#endif // PARALLEL

  data_ = BasicTimeAdvance(data_,sr_,rs_,eos_,dt,geometry_,
			   innerBC_, outerBC_);

  time_ += dt;
  cycle_++;
}

void SRHDSimulation::timeAdvance2ndOrder(void)
{
#ifdef PARALLEL
  throw("Not implemented yet");
#endif // PARALLEL
  const double dt = calcTimeStep();

  const NewHydroSnapshot<vector<double>, vector<Primitive> > mid = BasicTimeAdvance(data_,sr_,rs_,eos_,0.5*dt,
					     geometry_,
					     innerBC_,
					     outerBC_);

  CalcFluxes(mid,
	     sr_,rs_,dt,
	     innerBC_,
	     outerBC_,
	     psvs_);

  const vector<bool> filter = NeedUpdate(psvs_);

  update_new_conserved(psvs_,
		       data_.cells,
		       restMass_,
		       dt, 
		       geometry_,
		       data_.edges,
		       consVars_);

  UpdatePrimitives(consVars_, eos_, filter, data_.cells);

  time_ += dt;
  cycle_++;
}

namespace {

  double new_calc_time_step
  (const NewHydroSnapshot<vector<double>, vector<Primitive> >& data,
   const vector<RiemannSolution>& psvs,
   const EquationOfState& eos,
   const double cfl)
  {
    double its = 0; // Inverse time step
    for(size_t i=0;i<data.cells.size();++i){
      const double width = data.edges[i+1] - data.edges[i];
      const Primitive& cell = data.cells[i];
      const double ba = eos.dp2ba(cell.Density,
				  cell.Pressure);
      const double bc = celerity2velocity(data.cells[i].Celerity);
      const double bl = celerity2velocity(psvs[i].Celerity);
      const double br = celerity2velocity(psvs[i+1].Celerity);
      if(RelVelAdd(bc,ba)>br)
	its = fmax(its, (bc-br+ba*(1-bc*br))/(width*(1+bc*ba)));
      if(bl>RelVelAdd(bc,-ba))
	its = fmax(its, (bl-bc+ba*(1-bc*bl))/(width*(1-bc*ba)));
      its = fmax(its, 2*fabs(bl-br)/width);
    }
    return cfl/its;
  }
}

void SRHDSimulation::timeAdvance(void)
{
  //  double dt = TimeStep();

  CalcFluxes(data_,
	     sr_,rs_,0,
	     innerBC_,
	     outerBC_,
	     psvs_);

#ifdef PARALLEL
  const double dt_candidate =
    new_calc_time_step(data_,
		       psvs_,
		       eos_,
		       cfl_);
  double temp = dt_candidate;
  MPI_Allreduce(&dt_candidate,
		&temp,
		1,
		MPI_DOUBLE,
		MPI_MIN,
		MPI_COMM_WORLD);
  const double dt = temp;
  spdlog::debug("dt = {0}", dt);
#else
  const double dt = new_calc_time_step(data_,
				       psvs_,
				       eos_,
				       cfl_);
#endif // PARALLEL

  update_new_conserved(psvs_,
		       data_.cells,
		       restMass_,
		       dt, geometry_,
		       data_.edges,
		       consVars_);

  for(size_t i=1;i<data_.edges.size();++i)
    assert(data_.edges[i]>data_.edges[i-1]);
  for_each(consVars_.begin(),
	   consVars_.end(),
	   [](const NewConserved& cv)
	   {assert(cv.mass>0);});

  const vector<bool> filter = NeedUpdate(psvs_);

  UpdatePrimitives(consVars_, eos_, filter, data_.cells);

  time_ += dt;
  cycle_++;
}

void SRHDSimulation::calcConservedFromPrimitive(void)
{
  consVars_ = primitives_to_new_conserveds(data_.cells,eos_);
}

// Diagnostics

const decltype(SRHDSimulation::data_)&
SRHDSimulation::getHydroSnapshot(void) const
{
  return data_;
}

const vector<RiemannSolution>& 
SRHDSimulation::getRiemannSolutions(void) const
{
  return psvs_;
}

const vector<NewConserved>& SRHDSimulation::getConserved(void) const
{
  return consVars_;
}

double SRHDSimulation::getTime(void) const
{
  return time_;
}

int SRHDSimulation::getCycle(void) const
{
  return cycle_;
}

const EquationOfState& SRHDSimulation::getEOS(void) const
{
  return eos_;
}
