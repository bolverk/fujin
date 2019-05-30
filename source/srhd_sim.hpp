/*! \file srhd_sim.hpp
  \brief Relativistic, Godunov type, lagrangian 1d hydrodynamic simulation
  \author Almog Yalinewich
 */

#ifndef SRHD_SIM_HPP
#define SRHD_SIM_HPP 1

#include <vector>
#include "riemann_solver.hpp"
#include "spatial_distribution.hpp"
#include "equation_of_state.hpp"
#include "utilities.hpp"
#include "spatial_reconstruction.hpp"
#include "boundary_condition.hpp"
#include "geometry.hpp"

/*! \brief Special relativistic hydrodynamic simulation
  \details The simulation is based on <a href="http://adsabs.harvard.edu/abs/2000A%26A...358.1157D"> F. Daigne & R. Moskovitch, "Gamma-ray bursts from internal shocks in a relativistic wind: a hydrodynamical study", A&A, v. 358 p 1157-1166 (2000) </a>
 */
class SRHDSimulation
{
private:
  //! \brief Hydrodynamic state (edges and cells)
  HydroSnapshot data_;
  //! \brief Equation of state
  EquationOfState const& eos;
  //! \brief Riemann solutions
  vector<RiemannSolution> psvs;
  //! \brief rs Riemann sovler
  RiemannSolver const& rs;
  //! \brief Spatial reconstruction method
  SpatialReconstruction const& sr;
  //! \brief Courant Friedrichs Lewi factor
  double CourantFactor;
  //! \brief Vectors of conserved variables
  vector<NewConserved> ConsVars;
  //! \brief Rest masses of the cells
  vector<double> RestMass;
  //! \brief Geometry
  Geometry const& geometry_;
  //! \brief Virtual time
  double Time;
  //! \brief Cycle number
  int Cycle;
  //! \brief Innter boundary condition
  BoundaryCondition const& pInnerBC;
  //! \brief Outer boundary condition
  BoundaryCondition const& pOuterBC;

public:
  /*! \brief Class constructor
    \param vertices List of vertices
    \param density_distribution Initial density distribution
    \param pressure_distribution Initial pressure distribution
    \param proper_velocity_distribution Initial proper velocity distribution
    \param inner_bc Inner boundary conditions
    \param outer_bc Outer boundary conditions
    \param eos Equation of state
    \param riemann_solver Reference to riemann solver
    \param interpolation_method Pointer to spatial reconstruction
    \param geometry Geometry
   */
  SRHDSimulation(const vector<double>& vertices,
		 const SpatialDistribution& density_distribution,
		 const SpatialDistribution& pressure_distribution,
		 const SpatialDistribution& proper_velocity_distribution,
		 const BoundaryCondition& inner_bc,
		 const BoundaryCondition& outer_bc,
		 const EquationOfState& eos,
		 const RiemannSolver& riemann_solver,
		 SpatialReconstruction& interpolation_method,
		 const Geometry& geometry);

  //! \brief Advances the simulation in time
  void TimeAdvance1stOrder(void);

  //! \brief Advances the simulation in time
  void TimeAdvance2ndOrder(void);

  //! \brief Advances the simulation in time
  void TimeAdvance(void);

  /*! \brief Calculates the conserved variables from the primitives 
   */
  void CalcConservedFromPrimitive(void);

  /*! \brief override default Courant Friedrichs Lewy coefficient
    \details Determins the ratio between the actual time step and the maximal time step
    \param cfl_new New CFL coefficient
   */
  void OverrideCFL(double cfl_new);

  // Diagnostics

  /*! \brief Calculates the time step for a certain cell
    \param Index Cell index
    \return Time step
  */
  double TimeStepForCell(size_t Index) const;

  /*! \brief Return the volume bounded by a certain vertex
    \param Index Vertex index
    \return Volume
  */
  double GetVolume(size_t Index) const;

  /*! \brief Returns the hydrodynamic snapshot
    \return Hydrodynamic snapshot
   */
  const HydroSnapshot& getHydroSnapshot(void) const;

  /*! \brief Returns a list of Riemann solutions
    \return Riemann solutions
   */
  const vector<RiemannSolution>& getRiemannSolutions(void) const;

  /*! \brief Returns conserved variables struct
    \param i Cell index
    \return Conserved variables
   */
  NewConserved GetConserved(size_t i) const;

  /*! \brief returns the centre of a cell
    \param index Cell index
    \return Cell centre
   */
  double GetCellCentre(size_t index) const;

  /*! \brief Returns the rest mass of the cell
    \param Index Call index
    \return Rest mass
   */
  double GetCellRestMass(size_t Index) const;

  /*! \brief Calculates the time step
    \return Time step
  */
  double TimeStep(void) const;

  /*! \brief Calculates the area at a certain cell centre
    \param Index Cell index
    \return Area
   */
  double GetArea(size_t Index) const;

  /*! \brief Returns the time
  \return Time
  */
  double GetTime(void) const;

  /*! \brief Returns cycle number
    \return Cycle number
   */
  int GetCycle(void) const;

  /*! \brief Returns the rest mass of a cell
    \param index Cell index
    \return Rest mass
   */
  double GetRestMass(size_t index) const;

  /*! \brief Returns the equation of state
    \return Equation of state
   */
  const EquationOfState& getEOS(void) const;
};

#endif
