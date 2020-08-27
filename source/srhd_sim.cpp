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

double SRHDSimulation::GetCellCentre(size_t index) const
{
  return 0.5*(data_.edges[index]+data_.edges[index+1]);
}

double SRHDSimulation::GetArea(size_t Index) const
{
  return geometry_.calcArea(GetCellCentre(Index));
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
 SpatialDistribution const& DensityDistribution,
 SpatialDistribution const& PressureDistribution,
 SpatialDistribution const& VelocityDistribution,
 BoundaryCondition const& piInnerBC,
 BoundaryCondition const& piOuterBC,
 EquationOfState const& reos,
 RiemannSolver const& rRiemannSolver,
 SpatialReconstruction& pInterpolationMethod,
 Geometry const& geometry):
  data_(distribute_vertices(rVertices),
	InitCells(distribute_vertices(rVertices),
		  DensityDistribution,
		  PressureDistribution,
		  VelocityDistribution)),
  eos(reos), 
  psvs(distribute_vertices(rVertices).size(),
       RiemannSolution()),
  rs(rRiemannSolver), 
  sr(pInterpolationMethod),
  CourantFactor(1./3.), 
  ConsVars(primitives_to_new_conserveds(data_.cells,reos)),
  RestMass(serial_generate(RestMassCalculator(data_,geometry))),
  geometry_(geometry),
  Time(0),Cycle(0),
  pInnerBC(piInnerBC), pOuterBC(piOuterBC) {}

SRHDSimulation::SRHDSimulation
(const HydroSnapshot& init_cond,
 BoundaryCondition const& piInnerBC,
 BoundaryCondition const& piOuterBC,
 EquationOfState const& reos,
 RiemannSolver const& rRiemannSolver,
 SpatialReconstruction& pInterpolationMethod,
 Geometry const& geometry):
  data_(init_cond),
  eos(reos), 
  psvs(init_cond.edges.size(),
       RiemannSolution()),
  rs(rRiemannSolver), 
  sr(pInterpolationMethod),
  CourantFactor(1./3.), 
  ConsVars(primitives_to_new_conserveds(data_.cells,reos)),
  RestMass(serial_generate(RestMassCalculator(data_,geometry))),
  geometry_(geometry),
  Time(0),Cycle(0),
  pInnerBC(piInnerBC), pOuterBC(piOuterBC) {}

void SRHDSimulation::OverrideCFL(double cfl_new)
{
  CourantFactor = cfl_new;
}

double SRHDSimulation::TimeStep(void) const
{
  return CourantFactor*MaxTimeStep(data_.edges, data_.cells,eos);
}

double SRHDSimulation::TimeStepForCell(size_t i) const
{
  return MaxTimeStep(data_.edges[i+1]-data_.edges[i],
		     data_.cells[i],eos);
}

double SRHDSimulation::GetRestMass(size_t Index) const
{
  return RestMass[Index];
}

double SRHDSimulation::GetVolume(size_t Index) const
{
  return geometry_.calcVolume(data_.edges[Index]);
}

void SRHDSimulation::TimeAdvance1stOrder(void)
{
#ifdef PARALLEL
  const double dt_candidate = TimeStep();
  double temp = dt_candidate;
  MPI_Allreduce(&dt_candidate,&temp,1,
		MPI_DOUBLE,
		MPI_MIN,
		MPI_COMM_WORLD);
  const double dt = temp;
  spdlog::debug("dt = {0}",dt);
#else
  const double dt = TimeStep();
#endif // PARALLEL

  data_ = BasicTimeAdvance(data_,sr,rs,eos,dt,geometry_,
			   pInnerBC, pOuterBC);

  Time += dt;
  Cycle++;
}

void SRHDSimulation::TimeAdvance2ndOrder(void)
{
#ifdef PARALLEL
  throw("Not implemented yet");
#endif // PARALLEL
  const double dt = TimeStep();

  const HydroSnapshot mid = BasicTimeAdvance(data_,sr,rs,eos,0.5*dt,
					     geometry_,
					     pInnerBC,pOuterBC);

  CalcFluxes(mid,
	     sr,rs,dt,
	     pInnerBC,
	     pOuterBC,
	     psvs);

  const vector<bool> filter = NeedUpdate(psvs);

  /*
  UpdateConserved(psvs,
		  RestMass,
		  dt, 
		  geometry_,
		  data_.edges,
		  ConsVars);
  */
  update_new_conserved(psvs,
		       data_.cells,
		       RestMass,
		       dt, 
		       geometry_,
		       data_.edges,
		       ConsVars);

  UpdatePrimitives(ConsVars, eos, filter, data_.cells);

  Time += dt;
  Cycle++;
}

namespace {

  double new_calc_time_step(const HydroSnapshot& data,
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

void SRHDSimulation::TimeAdvance(void)
{
  //  double dt = TimeStep();

  CalcFluxes(data_,
	     sr,rs,0,
	     pInnerBC,
	     pOuterBC,
	     psvs);

#ifdef PARALLEL
  const double dt_candidate =
    new_calc_time_step(data_,
		       psvs,
		       eos,
		       CourantFactor);
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
				       psvs,
				       eos,
				       CourantFactor);
#endif // PARALLEL

  update_new_conserved(psvs,
		       data_.cells,
		       RestMass,
		       dt, geometry_,
		       data_.edges,
		       ConsVars);

  for(size_t i=1;i<data_.edges.size();++i)
    assert(data_.edges[i]>data_.edges[i-1]);
  for(size_t i=0;i<ConsVars.size();++i)
    assert(ConsVars[i].mass>0);

  const vector<bool> filter = NeedUpdate(psvs);

  UpdatePrimitives(ConsVars, eos, filter, data_.cells);

  Time += dt;
  Cycle++;
}

void SRHDSimulation::CalcConservedFromPrimitive(void)
{
  ConsVars = primitives_to_new_conserveds(data_.cells,eos);
}

// Diagnostics

const HydroSnapshot& SRHDSimulation::getHydroSnapshot(void) const
{
  return data_;
}

const vector<RiemannSolution>& 
SRHDSimulation::getRiemannSolutions(void) const
{
  return psvs;
}

NewConserved SRHDSimulation::GetConserved(size_t i) const
{
  return ConsVars[i];
}

double SRHDSimulation::GetCellRestMass(size_t i) const
{
  return RestMass[i];
}

double SRHDSimulation::GetTime(void) const
{
  return Time;
}

int SRHDSimulation::GetCycle(void) const
{
  return Cycle;
}

const EquationOfState& SRHDSimulation::getEOS(void) const
{
  return eos;
}
