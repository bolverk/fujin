#include <limits>
#include <cassert>
#include <algorithm>
#include "srhydro.hpp"
#include "utilities.hpp"
#include "advanced_hydrodynamic_variables.hpp"
#include "universal_error.hpp"
#ifdef PARALLEL
#include "mpi.h"
#endif // PARALLEL

using std::min_element;
using std::pair;

namespace {

  //! \brief Calculates averages of juxtaposed entries
  template<class T> class MidValues: public Index2Member<T>
  {
  public:

    /*! \brief Class constructor
      \param v Original array
     */
    explicit MidValues(const vector<T>& v):
      v_(v) {}

    size_t getLength(void) const 
    {
      return v_.size();
    }

    T operator()(size_t i) const
    {
      return 0.5*(v_[i]+v_[i+1]);
    }

  private:
    //! \brief Original array
    const vector<T>& v_;
  };

  /*! \brief Calculates the areas of cells
    \param vertices Position of the vertices
    \param geometry Geometry
    \return Areas of cells
  */
  vector<double> CellAreas(vector<double> const& vertices, 
			   Geometry const& geometry)
  {
    vector<double> res(vertices.size()-1);
    transform(vertices.begin(),
	      vertices.end()-1,
	      vertices.begin()+1,
	      res.begin(),
	      [&geometry](const double x, const double y)
	      {return geometry.calcArea(0.5*(x+y));});
    return res;
  }

  /*! \brief Volumes contained within vertices
    \param vertices Position of vertices
    \param geometry Geometry
    \return Volumes
  */
  vector<double> VerticesVolumes(vector<double> vertices, 
				 Geometry const& geometry)
  {
    vector<double> res(vertices.size());
    transform(vertices.begin(),
	      vertices.end(),
	      res.begin(),
	      [&geometry](const double r)
	      {return geometry.calcVolume(r);});
    return res;
  }
}

double CellVolume(double rl, double rr, Geometry const& geometry)
{
  return geometry.calcVolume(rr) - geometry.calcVolume(rl);
}

RestMassCalculator::RestMassCalculator
(const NewHydroSnapshot<simple_vector, simple_vector>& hs,
 const Geometry& geometry):
  hs_(hs), geometry_(geometry) {}

size_t RestMassCalculator::getLength(void) const
{
  return hs_.cells.size();
}

double RestMassCalculator::operator()(size_t i) const
{
  const double lf = celerity2lorentz_factor(hs_.cells[i].Celerity);
  const double vol = geometry_.calcVolume(hs_.edges[i+1])-
    geometry_.calcVolume(hs_.edges[i]);
  return hs_.cells[i].Density*vol*lf;
}

namespace {

  //! \brief Initialises the cells
  class CellGenerator: public Index2Member<Primitive>
  {
  public:

    /*! \brief Class constructor
      \param grid Computational grid
      \param density Density distribution
      \param pressure Pressure distribution
      \param celerity Celerity distribution
     */
    CellGenerator(const vector<double>& grid,
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
    const vector<double>& grid_;
    //! \brief Density distribution
    const SpatialDistribution& density_;
    //! \brief Pressure distribution
    const SpatialDistribution& pressure_;
    //! \brief Velocity distribution
    const SpatialDistribution& celerity_;
  };
}

vector<Primitive> InitCells(vector<double> const& v,
			    SpatialDistribution const& dd,
			    SpatialDistribution const& pd,
			    SpatialDistribution const& vd)
{
  return serial_generate(CellGenerator(v,dd,pd,vd));
}

namespace {

  //! \brief Converts primitives to conserved variables
  class Primitive2ConservedConverter: public Index2Member<Conserved>
  {
  public:

    /*! \brief Class constructor
      \param p_list Computational cells
      \param eos Equation of state
     */
    Primitive2ConservedConverter(const vector<Primitive>& p_list,
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
    const vector<Primitive>& p_list_;
    //! \brief Equation of state
    const EquationOfState& eos_;
  };

  //! \brief Converts primitives to conserved variables
  class Primitive2NewConservedConverter: public Index2Member<NewConserved>
  {
  public:

    /*! \brief Class constructor
      \param p_list Computational cells
      \param eos Equation of state
     */
    Primitive2NewConservedConverter(const vector<Primitive>& p_list,
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
    const vector<Primitive>& p_list_;
    //! \brief Equation of state
    const EquationOfState& eos_;
  };
}

vector<Conserved> Primitives2Conserveds
(vector<Primitive> const& p,
 const EquationOfState& eos )
{
  return serial_generate(Primitive2ConservedConverter(p,eos));
}

vector<NewConserved> primitives_to_new_conserveds
(const vector<Primitive>& p_list,
 const EquationOfState& eos)
{
  return serial_generate(Primitive2NewConservedConverter(p_list,eos));
}

double MaxTimeStep(double width, Primitive const& p, 
		   const EquationOfState& eos)
{
  const double sound_speed = eos.dp2ba(p.Density, p.Pressure);
  const double b_c = celerity2velocity(p.Celerity);
  const double g2_c = 1+pow(p.Celerity,2);
  return width*(1-sound_speed*fabs(b_c))*g2_c/sound_speed;
}

double MaxTimeStep(vector<double> const& v_list, 
		   vector<Primitive> const& p_list,
		   const EquationOfState& eos)
{
  //  return min(MaxTimeSteps(v,p,eos));
  //  return min(MaxTimeStepCalculator(v,p,eos));
  vector<double> valid_time_steps;
  for(size_t i=0;i<p_list.size();++i){
    const Primitive p = p_list[i];
    const double ba = eos.dp2ba(p.Density, p.Pressure);
    if(std::numeric_limits<double>::epsilon()<ba){
      const double width = v_list[i+1] - v_list[i];
      assert(width>std::numeric_limits<double>::epsilon());
      //      valid_time_steps.push_back(width/max_speed);
      const double bc = celerity2velocity(p.Celerity);
      const double g2c = 1.+pow(p.Celerity,2);
      valid_time_steps.push_back((1-fabs(bc)*ba)*width*g2c/ba);
    }
  }
  return *min_element(valid_time_steps.begin(), valid_time_steps.end());
}

namespace{
#ifdef PARALLEL
  vector<double> serialise_primitive
    (const Primitive& p)
  {
    vector<double> res(3,0);
    res.at(0) = p.Density;
    res.at(1) = p.Pressure;
    res.at(2) = p.Celerity;
    return res;
  }
  
  Primitive unserialise_primitive
    (const vector<double>& v)
  {
    assert(v.size()==3);
    Primitive res;
    res.Density = v.at(0);
    res.Pressure = v.at(1);
    res.Celerity = v.at(2);
    return res;
  }

  int get_mpi_rank(void)
  {
    int res;
    MPI_Comm_rank(MPI_COMM_WORLD, &res);
    return res;
  }

  int get_mpi_size(void)
  {
    int res;
    MPI_Comm_size(MPI_COMM_WORLD, &res);
    return res;
  }
#endif // PARALLEL
}

void CalcFluxes
(const NewHydroSnapshot<simple_vector, simple_vector>& data,
 const SpatialReconstruction& sr,
 const RiemannSolver& rs,
 double dt,
 const BoundaryCondition& lbc,
 const BoundaryCondition& rbc,
 vector<RiemannSolution>& psvs)
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
  psvs[0] = lbc.CalcRS(0,data.cells);
  psvs[psvs.size()-1] = rbc.CalcRS(data.edges.size()-1,data.cells);
#endif // PARALLEL
  const vector<std::pair<Primitive, Primitive> > interp_vals =
    sr.interpolateAll(data,dt);
  transform(interp_vals.begin(),
	    interp_vals.end(),
	    psvs.begin()+1,
	    [&rs](const pair<Primitive, Primitive>& pp)
	    {return rs(pp.first, pp.second);});
}

namespace{

  /*! \brief Calculates the fluxes
    \param data Hydrodynamic snapshot
    \param sr Spatial reconstruction
    \param rs Riemann solver
    \param dt Time step
    \param lbc Left boundary condition
    \param rbc Right boundary condition
    \return Fluxes
  */
  vector<RiemannSolution> CalcFluxes
    (const NewHydroSnapshot<simple_vector, simple_vector>& data,
     const SpatialReconstruction& sr,
     const RiemannSolver& rs,
     double dt,
     const BoundaryCondition& lbc,
     const BoundaryCondition& rbc)
  {
    vector<RiemannSolution> res(data.edges.size());
    vector<Primitive> new_cells(data.cells.size());
    CalcFluxes(data,sr,rs,dt,lbc,rbc,res);
    return res;
  }
}

void UpdateConserved(vector<RiemannSolution>  const& psvs,
		     vector<double> const& rest_mass,
		     double dt, 
		     Geometry const& geometry,
		     vector<double>& vertices,
		     vector<Conserved>& conserved)
{
  vector<double> volume_old = VerticesVolumes
    (vertices,geometry);
  vector<double> area_old = CellAreas(vertices,geometry);
  transform(psvs.begin(),
	    psvs.end(),
	    vertices.begin(),
	    vertices.begin(),
	    [&dt](const RiemannSolution& rsol, const double pos)
	    {return pos+dt*celerity2velocity(rsol.Celerity);});

  vector<double> volume_new = VerticesVolumes(vertices,geometry);
  vector<double> area_new = CellAreas(vertices,geometry);

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

namespace {

  /*! \brief Flux of the positive conserved variable
    \param p Pressure
    \param w Celerity
    \return Flux of the positive conserved variable
   */
  double positive_flux(double p, double w)
  {
    return p*(w>0 ?
	      1+celerity2velocity(w) :
	      1./(1+pow(w,2)-w*sqrt(1+pow(w,2))));
  }

  /*! \brief Calculates the flux of the positive conserved variable
    \param rs Riemann solution
    \return Flux of the positive conserved variable
   */
  double positive_flux(const RiemannSolution& rs)
  {
    return positive_flux(rs.Pressure, rs.Celerity);
  }

  /*! \brief Flux of the negative conserved variable
    \param p Pressure
    \param w Celerity
    \return Flux of the negative conserved variable
   */
  double negative_flux(double p, double w)
  {
    return p*(w<0 ?
	      -(1-celerity2velocity(w)) :
	      -1./(1+pow(w,2)+w*sqrt(1+pow(w,2))));
  }
     
  /*! \brief Calculates the flux of the negative conserved variable
    \param rs Riemann solution
    \return Flux of the negative conserved variable
   */ 
  double negative_flux(const RiemannSolution& rs)
  {
    return negative_flux(rs.Pressure, rs.Celerity);
  }
}

namespace{
  vector<double> calc_all_vertex_areas
    (const Geometry& geo,
     const vector<double>& vertices)
  {
    vector<double> res(vertices.size());
    transform(vertices.begin(),
	      vertices.end(),
	      res.begin(),
	      [&geo](const double r)
	      {return geo.calcArea(r);});
    return res;
  }
}

void update_new_conserved(const vector<RiemannSolution>& psvs,
			  const vector<Primitive>& cells,
			  const vector<double>& rest_mass,
			  double dt, 
			  const Geometry& geometry,
			  vector<double>& vertices,
			  vector<NewConserved>& conserved)
{
  const vector<double> vertex_areas = calc_all_vertex_areas(geometry, vertices);

  transform(psvs.begin(),
	    psvs.end(),
	    vertices.begin(),
	    vertices.begin(),
	    [&dt](const RiemannSolution& rsol, double pos)
	    {return pos + dt*celerity2velocity(rsol.Celerity);});

  const vector<double> volume_new = VerticesVolumes(vertices,geometry);

  for(size_t i=0;i<conserved.size();++i){
    conserved[i].positive += 
      (positive_flux(psvs[i])*vertex_areas[i]-
       positive_flux(psvs[i+1])*vertex_areas[i+1])*dt/
      rest_mass[i]+
      cells[i].Pressure*(vertex_areas[i+1]-vertex_areas[i])*dt/rest_mass[i];
    conserved[i].negative +=
      (negative_flux(psvs[i])*vertex_areas[i]-
       negative_flux(psvs[i+1])*vertex_areas[i+1])*dt/
      rest_mass[i]-
      cells[i].Pressure*(vertex_areas[i+1]-vertex_areas[i])*dt/rest_mass[i];
    conserved[i].mass = rest_mass[i]/
      (volume_new[i+1]-volume_new[i]);
  }
}

void UpdatePrimitives(vector<Conserved> const& conserved,
		      EquationOfState const& eos,
		      vector<Primitive>& cells)
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

void UpdatePrimitives(vector<NewConserved> const& conserved,
		      EquationOfState const& eos,
		      vector<Primitive>& cells)
{
  transform(conserved.begin(),
	    conserved.end(),
	    cells.begin(),
	    cells.begin(),
	    [&eos](const NewConserved& cons,
		   const Primitive& prim)
	    {return eos.Conserved2Primitive
		(prim, cons);});
}

vector<bool> NeedUpdate(vector<RiemannSolution> const& psvs)
{
  const class Checker: public Index2Member<bool>
  {
  public:
    explicit Checker(const vector<RiemannSolution>& psvs_i):
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
    const vector<RiemannSolution>& psvs_;
  } checker(psvs);
  return serial_generate(checker);
}

void UpdatePrimitives(vector<Conserved> const& conserved,
		      EquationOfState const& eos,
		      vector<bool> const& filter,
		      vector<Primitive>& cells)
{
  for(size_t i=0;i<cells.size();++i){
    if(filter[i])
      cells[i] = eos.Conserved2Primitive(cells[i],
					 old_to_new_conserved(conserved[i]));
  }
}

void UpdatePrimitives(vector<NewConserved> const& conserved,
		      EquationOfState const& eos,
		      vector<bool> const& filter,
		      vector<Primitive>& cells)
{
  for(size_t i=0;i<cells.size();++i){
    if(filter[i])
      cells[i] = eos.Conserved2Primitive(cells[i], conserved[i]);
  }
}

NewHydroSnapshot<simple_vector, simple_vector>
BasicTimeAdvance
(const NewHydroSnapshot<simple_vector, simple_vector>& data,
 SpatialReconstruction const& sr,
 RiemannSolver const& rs,
 EquationOfState const& eos,
 double dt, 
 Geometry const& geometry,
 BoundaryCondition const& lbc,
 BoundaryCondition const& rbc)
{
  NewHydroSnapshot<simple_vector, simple_vector> res(data);

  vector<RiemannSolution> fluxes=
    CalcFluxes(data,sr,rs,dt,lbc,rbc);

  const vector<double> rest_masses = serial_generate
    (RestMassCalculator(data, geometry));

  vector<Conserved> conserved = Primitives2Conserveds(data.cells,eos);
  UpdateConserved(fluxes,rest_masses,dt,geometry,
		  res.edges,conserved);

  vector<bool> filter = NeedUpdate(fluxes);
  UpdatePrimitives(conserved,eos,filter,res.cells);

  return NewHydroSnapshot<simple_vector, simple_vector>(res.edges, res.cells);
}

HydroSnapshot TimeAdvanceRK2(const HydroSnapshot& old,
			     const SpatialReconstruction& sr,
			     const RiemannSolver& rs,
			     const EquationOfState& eos,
			     double dt,
			     const Geometry& geometry,
			     const BoundaryCondition& lbc,
			     const BoundaryCondition& rbc)
{
  const NewHydroSnapshot<simple_vector, simple_vector> mid = BasicTimeAdvance
    (old, sr, rs, eos, 0.5*dt, geometry, lbc,rbc);
  
  const vector<RiemannSolution> fluxes = CalcFluxes
    (mid, sr, rs, dt, lbc, rbc);

  const vector<double> rest_masses = serial_generate
    (RestMassCalculator(old,geometry));
  vector<Conserved> conserved = Primitives2Conserveds(old.cells,eos);
  HydroSnapshot res = old;
  UpdateConserved(fluxes,rest_masses,dt,geometry,
		  res.edges,conserved);
  UpdatePrimitives(conserved,eos,res.cells);
  return res;
}
			     
