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

template<class T> using simple_vector = vector<T>;

/*! \brief Calculates the volume of the cell
  \param rl Radius of inner cell interface
  \param rr Radius of outer cell interface
  \param geometry Geometry
  \return Volume of cell
*/
double CellVolume(double rl, double rr,
		  Geometry const& geometry);

/*! \brief Calculates the rest masses
  \param v Vertices
  \param p Primitive variables
  \param geometry Geometry
  \return Vector of rest masses
*/
/*
vector<double> CalcRestMasses(const HydroSnapshot& hs,
			      const Geometry& geometry);
*/
class RestMassCalculator: public Index2Member<double>
{
public:

  /*! \brief Class constructor
    \param hs Hydrodynamic snapshot
    \param geometry Geometry
   */
  RestMassCalculator(const NewHydroSnapshot<simple_vector, simple_vector>& hs,
		     const Geometry& geometry);

  size_t getLength(void) const;

  double operator()(size_t i) const;
private:
  //! \brief Hydrodynamic snapshot
  const NewHydroSnapshot<simple_vector, simple_vector>& hs_;

  //! \brief Geometry
  const Geometry& geometry_;
};

/*! \brief Initialises primitive variables
  \param v Vertices
  \param dd Density distribution
  \param pd Pressure distribution
  \param vd Velocity distribution
  \return Vector of primitive variables
*/
vector<Primitive> InitCells(vector<double> const& v,
			    SpatialDistribution const& dd,
			    SpatialDistribution const& pd,
			    SpatialDistribution const& vd);

/*! \brief Calculates the vector of conserved variables
  \param p Primitive variables
  \param eos Equation of state
  \return Vector of conserved variables
*/
vector<Conserved> Primitives2Conserveds
(vector<Primitive> const& p, const EquationOfState& eos);

/*! \brief Calculates the vector of conserved variables
  \param p Primitive variables
  \param eos Equation of state
  \return Vector of conserved variables
*/
vector<NewConserved> primitives_to_new_conserveds
(vector<Primitive> const& p, const EquationOfState& eos);

/*! \brief Calculates the maximum time step according to CFL condition, for a single cell
  \param width Cell width
  \param p Primitive variables
  \param eos Equation of state
  \return Max time step
*/
double MaxTimeStepSingle(double width, Primitive const& p,
		   const EquationOfState& eos );

/*! \brief Calculates the maximum time step according to CFL condition, for the entire grid
  \param v Vertices
  \param p Cells
  \param eos Equation of state
  \return Max time step
*/
/*
double MaxTimeStep(vector<double> const& v, vector<Primitive> const& p,
		   const EquationOfState& eos);
*/

template<template<class> class CE, template<class> class CP>
double MaxTimeStep
(const CE<double>& v_list,
 const CP<Primitive>& p_list,
 const EquationOfState& eos)
{
  //  vector<double> valid_time_steps;
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
void CalcFluxes
(const NewHydroSnapshot<simple_vector, simple_vector>& data,
 const SpatialReconstruction& sr,
 const RiemannSolver& rs,
 double dt, 
 const BoundaryCondition& lbc,
 const BoundaryCondition& rbc,
 vector<RiemannSolution>& psvs);

/*! \brief Updates the conserved variables and vertices
  \param psvs Riemann solutions
  \param rest_mass Rest mass
  \param dt Time step
  \param geometry Geometry
  \param vertices Vertices
  \param conserved Conserved variables
*/
void UpdateConserved(vector<RiemannSolution>  const& psvs,
		     vector<double> const& rest_mass,
		     double dt, 
		     Geometry const& geometry,
		     vector<double>& vertices,
		     vector<Conserved>& conserved);

/*! \brief Updates the conserved variables and vertices
  \param psvs Riemann solutions
  \param cells Primitive variables
  \param rest_mass Rest mass
  \param dt Time step
  \param geometry Geometry
  \param vertices Vertices
  \param conserved Conserved variables
*/
void update_new_conserved(const vector<RiemannSolution>& psvs,
			  const vector<Primitive>& cells,
			  const vector<double>& rest_mass,
			  double dt, 
			  const Geometry& geometry,
			  vector<double>& vertices,
			  vector<NewConserved>& conserved);

/*! \brief Updates the primitive variables
  \param conserved Conserved variables
  \param eos Equation of state
  \param cells Primitive variables
*/
void UpdatePrimitives(vector<Conserved> const& conserved,
		      EquationOfState const& eos,
		      vector<Primitive>& cells);

/*! \brief Updates the primitive variables
  \param conserved Conserved variables
  \param eos Equation of state
  \param cells Primitive variables
*/
void UpdatePrimitives(vector<NewConserved> const& conserved,
		      EquationOfState const& eos,
		      vector<Primitive>& cells);

/*! \brief Updates the primitive variables
  \param conserved Conserved variables
  \param eos Equation of state
  \param filter array which determines which cells should be updated
  \param cells Primitive variables
*/
void UpdatePrimitives(vector<Conserved> const& conserved,
		      EquationOfState const& eos,
		      vector<bool> const& filter,
		      vector<Primitive>& cells);

/*! \brief Updates the primitive variables
  \param conserved Conserved variables
  \param eos Equation of state
  \param filter array which determines which cells should be updated
  \param cells Primitive variables
*/
void UpdatePrimitives(vector<NewConserved> const& conserved,
		      EquationOfState const& eos,
		      vector<bool> const& filter,
		      vector<Primitive>& cells);

/*! \brief Checks which cells need to be updated
  \details This can save run time by skipping cells with no velocity or pressure gradients
  \param psvs Riemann solutions at edges
  \return Array of boolean variables, denoting cells that should be updated by true
*/
vector<bool> NeedUpdate(vector<RiemannSolution> const& psvs);

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
NewHydroSnapshot<simple_vector, simple_vector> BasicTimeAdvance
(const NewHydroSnapshot<simple_vector, simple_vector>& data,
 const SpatialReconstruction& sr,
 const RiemannSolver& rs,
 const EquationOfState& eos,
 double dt, 
 const Geometry& geometry,
 const BoundaryCondition& lbc,
 const BoundaryCondition& rbc);

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
HydroSnapshot TimeAdvanceRK2(HydroSnapshot const& data,
			     SpatialReconstruction const& sr,
			     RiemannSolver const& rs,
			     EquationOfState const& eos,
			     double dt, 
			     Geometry const& geometry,
			     BoundaryCondition const& lbc,
			     BoundaryCondition const& rbc);

#endif // SRHYDRO_HPP
