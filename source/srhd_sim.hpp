/*! \file srhd_sim.hpp
  \brief Relativistic, Godunov type, lagrangian 1d hydrodynamic simulation
  \author Almog Yalinewich
*/

#ifndef SRHD_SIM_HPP
#define SRHD_SIM_HPP 1

#define SCAFFOLDING 0

#include <vector>
#include "riemann_solver.hpp"
#include "spatial_distribution.hpp"
#include "equation_of_state.hpp"
#include "utilities.hpp"
#include "spatial_reconstruction.hpp"
#include "boundary_condition.hpp"
#include "geometry.hpp"
#include <cassert>
#include "srhydro.hpp"

template<class CE, class CP> double new_calc_time_step
(const NewHydroSnapshot<CE, CP>& data,
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

vector<double> distribute_vertices1
(const vector<double>& vertices);

/*! \brief Special relativistic hydrodynamic simulation
  \details The simulation is based on <a href="http://adsabs.harvard.edu/abs/2000A%26A...358.1157D"> F. Daigne & R. Moskovitch, "Gamma-ray bursts from internal shocks in a relativistic wind: a hydrodynamical study", A&A, v. 358 p 1157-1166 (2000) </a>
*/
template<class CE, class CP> class SRHDSimulation
  {
  private:
    //! \brief Hydrodynamic state (edges and cells)
    NewHydroSnapshot<vector<double>, vector<Primitive> > data_;
    //! \brief Equation of state
    const EquationOfState& eos_;
    //! \brief Riemann solutions
    vector<RiemannSolution> psvs_;
    //! \brief rs Riemann sovler
    const RiemannSolver& rs_;
    //! \brief Spatial reconstruction method
    const SpatialReconstruction& sr_;
    //! \brief Courant Friedrichs Lewy factor
    double cfl_;
    //! \brief Vectors of conserved variables
    vector<NewConserved> consVars_;
    //! \brief Rest masses of the cells
    vector<double> restMass_;
    //! \brief Geometry
    const Geometry& geometry_;
    //! \brief Virtual time
    double time_;
    //! \brief Cycle number
    int cycle_;
    //! \brief Innter boundary condition
    const BoundaryCondition& innerBC_;
    //! \brief Outer boundary condition
    const BoundaryCondition& outerBC_;

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
		   const SpatialReconstruction& interpolation_method,
		   const Geometry& geometry):
    data_(distribute_vertices1(vertices),
	  InitCells(distribute_vertices1(vertices),
		    density_distribution,
		    pressure_distribution,
		    proper_velocity_distribution)),
      eos_(eos), 
      psvs_(distribute_vertices1(vertices).size(),
	    RiemannSolution()),
      rs_(riemann_solver), 
      sr_(interpolation_method),
      cfl_(1./3.), 
      consVars_(primitives_to_new_conserveds(data_.cells,eos)),
      restMass_(serial_generate(RestMassCalculator(data_,geometry))),
      geometry_(geometry),
      time_(0),
      cycle_(0),
      innerBC_(inner_bc), 
      outerBC_(outer_bc) {}

    /*! \brief Class constructor
      \param init_cond Initial conditions
      \param inner_bc Inner boundary conditions
      \param outer_bc Outer boundary conditions
      \param eos Equation of state
      \param rs Riemann solver
      \param interpolation_method Pointer to spatial reconstruction
      \param geometry Geometry
    */
  SRHDSimulation
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
  
    //! \brief Advances the simulation in time
  void timeAdvance1stOrder(void)
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

    //! \brief Advances the simulation in time
  void timeAdvance2ndOrder(void)
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

    //! \brief Advances the simulation in time
    void timeAdvance(void)
    {
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

    /*! \brief Calculates the conserved variables from the primitives 
     */
    void calcConservedFromPrimitive(void);

    /*! \brief override default Courant Friedrichs Lewy coefficient
      \details Determins the ratio between the actual time step and the maximal time step
      \param cfl_new New CFL coefficient
    */
    void overrideCFL(double cfl_new);

    // Diagnostics

    /*! \brief Calculates the time step for a certain cell
      \param Index Cell index
      \return Time step
    */
  double calcTimeStepForCell(size_t i)
  {
    return cfl_*MaxTimeStep(data_.edges, data_.cells,eos_);
  }

    /*! \brief Return the volume bounded by a certain vertex
      \param Index Vertex index
      \return Volume
    */
  double getVolume(size_t i) const
  {
    return geometry_.calcVolume(data_.edges.at(i));
  }

    /*! \brief Returns the hydrodynamic snapshot
      \return Hydrodynamic snapshot
    */
    const decltype(data_)& getHydroSnapshot(void) const
    {
      return data_;
    }

    /*! \brief Returns a list of Riemann solutions
      \return Riemann solutions
    */
    const vector<RiemannSolution>& getRiemannSolutions(void) const;

    /*! \brief Returns the conserved variables
      \return Conserved variables
    */
  const vector<NewConserved>& getConserved(void) const
  {
    return consVars_;
  }

    /*! \brief Calculates the time step
      \return Time step
    */
  double calcTimeStep(void) const
  {
    return cfl_*MaxTimeStep(data_.edges, data_.cells,eos_);
  }

    /*! \brief Calculates the area at a certain cell centre
      \param Index Cell index
      \return Area
    */
  double getArea(size_t i) const
  {
    const double r = 0.5*(data_.edges.at(i)+data_.edges.at(i+1));
    return geometry_.calcArea(r);
  }

    /*! \brief Returns the time
      \return Time
    */
    double getTime(void) const
    {
      return time_;
    }

    /*! \brief Returns cycle number
      \return Cycle number
    */
    int getCycle(void) const
    {
      return cycle_;
    }

    /*! \brief Returns the rest masses
      \return Rest masses
    */
  const vector<double>&  getRestMasses(void) const
  {
    return restMass_;
  }

    /*! \brief Returns the equation of state
      \return Equation of state
    */
    const EquationOfState& getEOS(void) const;
  };

#endif
