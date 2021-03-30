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
}

double CellVolume(double rl, double rr, Geometry const& geometry)
{
  return geometry.calcVolume(rr) - geometry.calcVolume(rl);
}

double MaxTimeStepSingle(double width, Primitive const& p, 
		   const EquationOfState& eos)
{
  const double sound_speed = eos.dp2ba(p.Density, p.Pressure);
  const double b_c = celerity2velocity(p.Celerity);
  const double g2_c = 1+pow(p.Celerity,2);
  return width*(1-sound_speed*fabs(b_c))*g2_c/sound_speed;
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

/*! \brief Flux of the positive conserved variable
  \param p Pressure
  \param w Celerity
  \return Flux of the positive conserved variable
*/
double srhydro::positive_flux(double p, double w)
{
  return p*(w>0 ?
	    1+celerity2velocity(w) :
	    1./(1+pow(w,2)-w*sqrt(1+pow(w,2))));
}

/*! \brief Calculates the flux of the positive conserved variable
  \param rs Riemann solution
  \return Flux of the positive conserved variable
*/
double srhydro::positive_flux(const RiemannSolution& rs)
{
  return srhydro::positive_flux(rs.Pressure, rs.Celerity);
}

double srhydro::negative_flux(const RiemannSolution& rs)
{
  return negative_flux(rs.Pressure, rs.Celerity);
}

  /*! \brief Flux of the negative conserved variable
    \param p Pressure
    \param w Celerity
    \return Flux of the negative conserved variable
   */
double srhydro::negative_flux(double p, double w)
{
  return p*(w<0 ?
	    -(1-celerity2velocity(w)) :
	    -1./(1+pow(w,2)+w*sqrt(1+pow(w,2))));
}
