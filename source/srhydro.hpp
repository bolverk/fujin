/*! \file srhydro.hpp
  \brief Hydrodynamic functions
  \author Almog Yalinewich
*/

#ifndef SRHYDRO_HPP
#define SRHYDRO_HPP 1

#include <vector>
#include "hydrodynamic_variables.hpp"
#include "spatial_distribution.hpp"
#include "equation_of_state.hpp"
#include "spatial_reconstruction.hpp"
#include "riemann_solver.hpp"
#include "boundary_condition.hpp"
#include "geometry.hpp"
#include "utilities.hpp"
#include <cassert>
#include "advanced_hydrodynamic_variables.hpp"

/*! \brief Calculates the volume of the cell
  \param rl Radius of inner cell interface
  \param rr Radius of outer cell interface
  \param geometry Geometry
  \return Volume of cell
*/
double CellVolume(double rl, double rr,
		  Geometry const& geometry);

template<template<class> class CE, template<class> class CP>
class RestMassCalculator: public Index2Member<double>
{
public:

  /*! \brief Class constructor
    \param hs Hydrodynamic snapshot
    \param geometry Geometry
  */
  RestMassCalculator(const NewHydroSnapshot<CE, CP>& hs,
		     const Geometry& geometry):
    hs_(hs), geometry_(geometry) {}

  size_t getLength(void) const
  {
    return hs_.cells.size();
  }

  double operator()(size_t i) const
  {
    const double lf = celerity2lorentz_factor(hs_.cells[i].Celerity);
    const double vol = geometry_.calcVolume(hs_.edges[i+1])-
      geometry_.calcVolume(hs_.edges[i]);
    return hs_.cells[i].Density*vol*lf;
  }
private:
  //! \brief Hydrodynamic snapshot
  const NewHydroSnapshot<CE, CP>& hs_;

  //! \brief Geometry
  const Geometry& geometry_;
};

namespace srhydro{
  //! \brief Initialises the cells
  template<template<class> class CE>
  class CellGenerator: public Index2Member<Primitive>
  {
  public:

    /*! \brief Class constructor
      \param grid Computational grid
      \param density Density distribution
      \param pressure Pressure distribution
      \param celerity Celerity distribution
     */
    CellGenerator(const CE<double>& grid,
		  const SpatialDistribution& density,
		  const SpatialDistribution& pressure,
		  const SpatialDistribution& celerity):
      grid_(grid),
      density_(density),
      pressure_(pressure),
      celerity_(celerity) {}

    size_t getLength(void) const
    {
      return grid_.size()-1;
    }
    
    Primitive operator()(size_t i) const
    {
      const double x = 0.5*(grid_[i]+grid_[i+1]);
      return Primitive(density_(x),
		       pressure_(x),
		       celerity_(x));
    }

  private:
    //! \brief Computational grid
    const CE<double>& grid_;
    //! \brief Density distribution
    const SpatialDistribution& density_;
    //! \brief Pressure distribution
    const SpatialDistribution& pressure_;
    //! \brief Velocity distribution
    const SpatialDistribution& celerity_;
  };
}

/*! \brief Initialises primitive variables
  \param v Vertices
  \param dd Density distribution
  \param pd Pressure distribution
  \param vd Velocity distribution
  \return List of primitive variables
*/
template<template<class> class CE, template<class> class CP>
CP<Primitive> InitCells(const CE<double>& v,
			const SpatialDistribution& dd,
			const SpatialDistribution& pd,
			const SpatialDistribution& vd)
{
  return serial_generate<Primitive, CP>(srhydro::CellGenerator<CE>(v,dd,pd,vd));
}

namespace srhydro
{
  template<template<class> class CP>
  class Primitive2NewConservedConverter: public Index2Member<NewConserved>
  {
  public:

    /*! \brief Class constructor
      \param p_list Computational cells
      \param eos Equation of state
    */
    Primitive2NewConservedConverter
    (const CP<Primitive>& p_list,
     const EquationOfState& eos):
      p_list_(p_list), eos_(eos) {}

    size_t getLength(void) const
    {
      return p_list_.size();
    }

    NewConserved operator()(size_t i) const
    {
      return primitive_to_new_conserved(p_list_[i],eos_);
    }

  private:
    //! \brief Computational cells
    const CP<Primitive>& p_list_;
    //! \brief Equation of state
    const EquationOfState& eos_;
  };
}

/*! \brief Calculates the list of conserved variables
  \param p Primitive variables
  \param eos Equation of state
  \return List of conserved variables
*/
template<template<class> class CP>
CP<NewConserved> primitives_to_new_conserveds
(const CP<Primitive>& p, const EquationOfState& eos)
{
  return serial_generate<NewConserved, CP>
    (srhydro::Primitive2NewConservedConverter<CP>(p,eos));
}

/*! \brief Calculates the maximum time step according to CFL condition, for a single cell
  \param width Cell width
  \param p Primitive variables
  \param eos Equation of state
  \return Max time step
*/
double MaxTimeStepSingle(double width, Primitive const& p,
			 const EquationOfState& eos );

template<template<class> class CE, template<class> class CP>
double MaxTimeStep
(const CE<double>& v_list,
 const CP<Primitive>& p_list,
 const EquationOfState& eos)
{
  double res = 0;
  for(size_t i=0;i<p_list.size();++i){
    const Primitive p = p_list[i];
    const double ba = eos.dp2ba(p.Density, p.Pressure);
    if(std::numeric_limits<double>::epsilon()<ba){
      const double width = v_list[i+1] - v_list[i];
      assert(width>std::numeric_limits<double>::epsilon());
      //      valid_time_steps.push_back(width/max_speed);
      const double bc = celerity2velocity(p.Celerity);
      const double g2c = 1.+pow(p.Celerity,2);
      //      valid_time_steps.push_back((1-fabs(bc)*ba)*width*g2c/ba);
      res = fmax(res, 1/((1-fabs(bc)*ba)*width*g2c/ba));
    }
  }
  //  return *min_element(valid_time_steps.begin(), valid_time_steps.end());
  return 1.0/res;
}

/*! \brief Calculates the hydrodyanmic fluxes
  \param data Hydrodynamic data (grid + cells)
  \param sr Spatial reconstruction
  \param rs Riemann solver
  \param dt Time step
  \param lbc Left boundary conditions
  \param rbc Right boundary conditions
  \param psvs Riemann solutions
*/
template<template<class> class CE, template<class> class CP> void CalcFluxes
(const NewHydroSnapshot<CE, CP>& data,
 const SpatialReconstruction<CE, CP>& sr,
 const RiemannSolver& rs,
 double dt, 
 const BoundaryCondition<CP>& lbc,
 const BoundaryCondition<CP>& rbc,
 CE<RiemannSolution>& psvs)
{
#ifdef PARALLEL
  MPI_Request send_left, send_right;
  // Send data
  if(get_mpi_rank()>0)
    MPI_Isend(&serialise_primitive(data.cells.front()).front(),
	      3,
	      MPI_DOUBLE,
	      get_mpi_rank()-1,
	      0,
	      MPI_COMM_WORLD,
	      &send_left);
  if(get_mpi_rank()<get_mpi_size()-1)
    MPI_Isend(&serialise_primitive(data.cells.back()).front(),
	      3,
	      MPI_DOUBLE,
	      get_mpi_rank()+1,
	      1,
	      MPI_COMM_WORLD,
	      &send_right);

  // Receive data and calculate edge fluxes
  if(get_mpi_rank()==0)
    psvs.front() = lbc.CalcRS(0,data.cells);
  else{
    vector<double> buffer(3,0);
    MPI_Recv(&buffer.front(),
	     3,
	     MPI_DOUBLE,
	     get_mpi_rank()-1,
	     1,
	     MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);
    MPI_Wait(&send_left, MPI_STATUS_IGNORE);
    psvs.front() = rs(unserialise_primitive(buffer),
		      data.cells.front());
  } 
  if(get_mpi_rank()==get_mpi_size()-1)
    psvs.back() = rbc.CalcRS(data.edges.size()-1,data.cells);
  else{
    vector<double> buffer(3,0);
    MPI_Recv(&buffer.front(),
	     3,
	     MPI_DOUBLE,
	     get_mpi_rank()+1,
	     0,
	     MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);
    MPI_Wait(&send_right, MPI_STATUS_IGNORE);
    psvs.back() = rs(data.cells.back(),
		     unserialise_primitive(buffer));
    
  }
#else
  psvs.front() = lbc(false,data.cells);
  psvs.back() = rbc(true,data.cells);
#endif // PARALLEL
  const CE<std::pair<Primitive, Primitive> > interp_vals =
    sr.interpolateAll(data,dt);
  transform(interp_vals.begin()+1,
	    interp_vals.end()-1,
	    psvs.begin()+1,
	    [&rs](const pair<Primitive, Primitive>& pp)
	    {return rs(pp.first, pp.second);});
}

namespace srhydro
{
  template<template<class> class CE, template<class> class CP>
  CE<RiemannSolution> CalcFluxes1
  (const NewHydroSnapshot<CE, CP>& data,
   const SpatialReconstruction<CE, CP>& sr,
   const RiemannSolver& rs,
   double dt,
   const BoundaryCondition<CP>& lbc,
   const BoundaryCondition<CP>& rbc)
  {
    CE<RiemannSolution> res;
    resize_if_necessary(res, data.edges.size());
    CalcFluxes<CE, CP>(data,sr,rs,dt,lbc,rbc,res);
    return res;
  }

  template<template<class> class CE>
  CE<double> VerticesVolumes(const CE<double>& vertices, 
			     const Geometry& geometry)
  {
    CE<double> res;
    resize_if_necessary(res, vertices.size());
    transform(vertices.begin(),
	      vertices.end(),
	      res.begin(),
	      [&geometry](const double r)
	      {return geometry.calcVolume(r);});
    return res;
  }

  template<template<class> class CE, template<class> class CP>
  CP<double> CellAreas
  (CE<double> const& vertices, 
   Geometry const& geometry)
  {
    CP<double> res;
    resize_if_necessary(res, vertices.size()-1);
    transform(vertices.begin(),
	      vertices.end()-1,
	      vertices.begin()+1,
	      res.begin(),
	      [&geometry](const double x, const double y)
	      {return geometry.calcArea(0.5*(x+y));});
    return res;
  }
}

/*! \brief Updates the conserved variables and vertices
  \param psvs Riemann solutions
  \param rest_mass Rest mass
  \param dt Time step
  \param geometry Geometry
  \param vertices Vertices
  \param conserved Conserved variables
*/
template<template<class> class CE, template<class> class CP>
void UpdateConserved(CE<RiemannSolution>  const& psvs,
		     CP<double> const& rest_mass,
		     double dt, 
		     Geometry const& geometry,
		     CE<double>& vertices,
		     CP<Conserved>& conserved)
{
  const CE<double> volume_old = srhydro::VerticesVolumes<CE>
    (vertices,geometry);
  const CE<double> area_old = srhydro::CellAreas<CE, CP>(vertices,geometry);
  transform(psvs.begin(),
	    psvs.end(),
	    vertices.begin(),
	    vertices.begin(),
	    [&dt](const RiemannSolution& rsol, const double pos)
	    {return pos+dt*celerity2velocity(rsol.Celerity);});

  const CE<double> volume_new = srhydro::VerticesVolumes<CE>(vertices,geometry);
  const CP<double> area_new = srhydro::CellAreas<CE, CP>(vertices,geometry);

  for(size_t i=0;i<conserved.size();++i){
    conserved[i].Energy += 
      (psvs[i].Pressure*(volume_new[i]-volume_old[i])-
       psvs[i+1].Pressure*(volume_new[i+1]-volume_old[i+1]))/
      rest_mass[i];
    conserved[i].Momentum +=
      (psvs[i].Pressure-psvs[i+1].Pressure)*dt*
      0.5*(area_new[i]+area_old[i])/rest_mass[i];
    conserved[i].Mass = rest_mass[i]/
      (volume_new[i+1]-volume_new[i]);
  }
}

namespace srhydro{

  double positive_flux(double p, double w);
      
  double positive_flux(const RiemannSolution& rs);

  double negative_flux(double p, double w);

  double negative_flux(const RiemannSolution& rs);

  template<template<class> class CE>
  CE<double> calc_all_vertex_areas
  (const Geometry& geo,
   const CE<double>& vertices)
  {
    CE<double> res;
    resize_if_necessary(res, vertices.size());
    transform(vertices.begin(),
	      vertices.end(),
	      res.begin(),
	      [&geo](const double r)
	      {return geo.calcArea(r);});
    return res;
  }

  template<template<class> class CP>
  class Primitive2ConservedConverter: public Index2Member<Conserved>
  {
  public:

    /*! \brief Class constructor
      \param p_list Computational cells
      \param eos Equation of state
     */
    Primitive2ConservedConverter
    (const CP<Primitive>& p_list,
     const EquationOfState& eos):
      p_list_(p_list), eos_(eos) {}

    size_t getLength(void) const
    {
      return p_list_.size();
    }

    Conserved operator()(size_t i) const
    {
      return Primitive2Conserved(p_list_[i],eos_);
    }

  private:
    //! \brief Computational cells
    const CP<Primitive>& p_list_;
    //! \brief Equation of state
    const EquationOfState& eos_;
  };

  template<template<class> class CP>
  CP<Conserved> Primitives2Conserveds
  (const CP<Primitive>& p,
   const EquationOfState& eos)
  {
    return serial_generate<Conserved, CP>
      (Primitive2ConservedConverter<CP>(p,eos));
  }
}

/*! \brief Updates the conserved variables and vertices
  \param psvs Riemann solutions
  \param cells Primitive variables
  \param rest_mass Rest mass
  \param dt Time step
  \param geometry Geometry
  \param vertices Vertices
  \param conserved Conserved variables
*/
template<template<class> class CE, template<class> class CP>
void update_new_conserved(const CE<RiemannSolution>& psvs,
			  const CP<Primitive>& cells,
			  const CP<double>& rest_mass,
			  double dt, 
			  const Geometry& geometry,
			  CE<double>& vertices,
			  CP<NewConserved>& conserved)
{
  const CE<double> vertex_areas =
    srhydro::calc_all_vertex_areas<CE>(geometry, vertices);
  
  transform(psvs.begin(),
	    psvs.end(),
	    vertices.begin(),
	    vertices.begin(),
	    [&dt](const RiemannSolution& rsol, double pos)
	    {return pos + dt*celerity2velocity(rsol.Celerity);});

  const CE<double> volume_new = srhydro::VerticesVolumes<CE>(vertices,geometry);

  for(size_t i=0;i<conserved.size();++i){
    conserved[i].positive += 
      (srhydro::positive_flux(psvs[i])*vertex_areas[i]-
       srhydro::positive_flux(psvs[i+1])*vertex_areas[i+1])*dt/
      rest_mass[i]+
      cells[i].Pressure*(vertex_areas[i+1]-vertex_areas[i])*dt/rest_mass[i];
    conserved[i].negative +=
      (srhydro::negative_flux(psvs[i])*vertex_areas[i]-
       srhydro::negative_flux(psvs[i+1])*vertex_areas[i+1])*dt/
      rest_mass[i]-
      cells[i].Pressure*(vertex_areas[i+1]-vertex_areas[i])*dt/rest_mass[i];
    conserved[i].mass = rest_mass[i]/
      (volume_new[i+1]-volume_new[i]);
  }
}

/*! \brief Updates the primitive variables
  \param conserved Conserved variables
  \param eos Equation of state
  \param cells Primitive variables
*/
template<template<class> class CP>
void UpdatePrimitives
(const CP<Conserved>& conserved,
 const EquationOfState& eos,
 CP<Primitive>& cells)
{
    transform(conserved.begin(),
	    conserved.end(),
	    cells.begin(),
	    cells.begin(),
	    [&eos](const Conserved& cons,
	       const Primitive& prim)
	    {return eos.Conserved2Primitive
		(prim,
		 old_to_new_conserved(cons));});
}

/*! \brief Updates the primitive variables
  \param conserved Conserved variables
  \param eos Equation of state
  \param filter array which determines which cells should be updated
  \param cells Primitive variables
*/
template<template<class> class CP>
void UpdatePrimitives(const CP<NewConserved>& conserved,
		      const EquationOfState& eos,
		      const CP<bool>& filter,
		      CP<Primitive>& cells)
{
  for(size_t i=0;i<cells.size();++i){
    if(filter[i])
      cells[i] = eos.Conserved2Primitive(cells[i], conserved[i]);
  }
}

template<template<class> class CP>
void UpdatePrimitives(CP<Conserved> const& conserved,
		      EquationOfState const& eos,
		      CP<bool> const& filter,
		      CP<Primitive>& cells)
{
  for(size_t i=0;i<cells.size();++i){
    if(filter[i])
      cells[i] = eos.Conserved2Primitive(cells[i],
					 old_to_new_conserved(conserved[i]));
  }
}

/*! \brief Checks which cells need to be updated
  \details This can save run time by skipping cells with no velocity or pressure gradients
  \param psvs Riemann solutions at edges
  \return Array of boolean variables, denoting cells that should be updated by true
*/
template<template<class> class CE, template<class> class CP>
CP<bool> NeedUpdate(CE<RiemannSolution> const& psvs)
{
  const class Checker: public Index2Member<bool>
  {
  public:
    explicit Checker(const CE<RiemannSolution>& psvs_i):
      psvs_(psvs_i) {}

    size_t getLength(void) const
    {
      return psvs_.size()-1;
    }

    bool operator()(size_t i) const
    {
      return !(effectively_zero(psvs_[i].Celerity)&&
	       effectively_zero(psvs_[i+1].Celerity)&&
	       effectively_zero(psvs_[i].Pressure-
				psvs_[i+1].Pressure));
    }

  private:
    const CE<RiemannSolution>& psvs_;
  } checker(psvs);
  return serial_generate<bool, CP>(checker);
}

/*! \brief Returns the hydrodynamic data in the next time step
  \param data Old hydrodynamic data
  \param sr Spatial reconstruction
  \param rs Riemann solver
  \param eos Equation of state
  \param dt Time step
  \param geometry Geometry of the problem
  \param lbc Left boundary conditions
  \param rbc Right boundary conditions
  \return Hydrodynamic data at the end of the time step
*/
template<template<class> class CE, template<class> class CP>
NewHydroSnapshot<CE, CP> BasicTimeAdvance
(const NewHydroSnapshot<CE, CP>& data,
 const SpatialReconstruction<CE, CP>& sr,
 const RiemannSolver& rs,
 const EquationOfState& eos,
 double dt, 
 const Geometry& geometry,
 const BoundaryCondition<CP>& lbc,
 const BoundaryCondition<CP>& rbc)
{
  NewHydroSnapshot<CE, CP> res(data);

  CE<RiemannSolution> fluxes=
    srhydro::CalcFluxes1<CE, CP>(data,sr,rs,dt,lbc,rbc);

  const CP<double> rest_masses = serial_generate
    (RestMassCalculator<CE, CP>(data, geometry));

  CP<Conserved> conserved = srhydro::Primitives2Conserveds<CP>(data.cells,eos);
  UpdateConserved<CE, CP>(fluxes,rest_masses,dt,geometry,
		  res.edges,conserved);

  const CP<bool> filter = NeedUpdate<CE, CP>(fluxes);
  UpdatePrimitives<CP>(conserved,eos,filter,res.cells);

  return NewHydroSnapshot<CE, CP>(res.edges, res.cells);
}

/*! \brief Second order time advance
  \param data Hydrodynamic snapshot
  \param sr Interpolation
  \param rs Riemann solver
  \param eos Equation of state
  \param dt Time step
  \param geometry Geometry
  \param lbc Left boundary condition
  \param rbc Right boundary condition
  \return Hydrodynamic snapshot
*/
/*
  HydroSnapshot TimeAdvanceRK2
  (HydroSnapshot const& data,
  SpatialReconstruction<spatialreconstruction const& sr,
  RiemannSolver const& rs,
  EquationOfState const& eos,
  double dt, 
  Geometry const& geometry,
  BoundaryCondition const& lbc,
  BoundaryCondition const& rbc);
*/

#endif // SRHYDRO_HPP
