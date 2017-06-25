#include <cmath>
#include "van_leer.hpp"
#include "universal_error.hpp"
#include "utilities.hpp"
#include "advanced_hydrodynamic_variables.hpp"

namespace {

  /*! \brief returns the sign of a number
    \param x Number
    \return Sign
  */
  double sgn(double x)
  {
    return (x>0) - (x<0);
  }

  /*! \brief Minmod
    \param x First argument
    \param y Second argument
    \return Minmod of x and y
  */
  double minmod(double x, double y)
  {
    return 0.5*(sgn(x)+sgn(y))*
      min(fabs(x),fabs(y));
  }

  /*! \brief Applies minmod filter to all members
    \param p1 First set of primitive variables
    \param p2 Second set of primitive variables
    \return Minmod filtered primitve variables 
   */
  Primitive minmod(Primitive const& p1,
		   Primitive const& p2)
  {
    return Primitive(minmod(p1.Density,p2.Density),
		     minmod(p1.Pressure,p2.Pressure),
		     minmod(p1.Celerity,p2.Celerity));
  }

  /*! \brief Calculates the slopes of the cells
    \param hs Hydrodynamic snapshot
    \return Derivatives of the hydrodynamic variables on the edges
   */
  vector<Primitive> calc_edge_slopes(const HydroSnapshot& hs)
  {
    const class DerivativeInterface: public Index2Member<Primitive>
    {
    public:

      DerivativeInterface(const HydroSnapshot& hs_i):
	hs_(hs_i) {}

      size_t getLength(void) const
      {
	return hs_.edges.size();
      }

      Primitive operator()(size_t i) const
      {
	if(i==0||i==hs_.edges.size()-1)
	  return Primitive();
	else
	  return (hs_.cells[i]-hs_.cells[i-1])/
	    (0.5*(hs_.edges[i+1]-hs_.edges[i-1]));
      }

    private:
      const HydroSnapshot& hs_;
    } deriv_interface(hs);
    return serial_generate(deriv_interface);
  }

  /*! \brief Calculates the slopes on the computational cells
    \param hs Hydrodynamic snapshot
    \return List of slopes
   */
  vector<Primitive> calc_cell_slopes(const HydroSnapshot& hs)
  {
    const vector<Primitive> edge_slopes = calc_edge_slopes(hs);
    const class Minmoder: public Index2Member<Primitive>
    {
    public:

      Minmoder(const vector<Primitive>& edge_slopes_i):
	edge_slopes_(edge_slopes_i) {}

      size_t getLength(void) const
      {
	return edge_slopes_.size()-1;
      }

      Primitive operator()(size_t i) const
      {
	return minmod(edge_slopes_[i],
		      edge_slopes_[i+1]);
      }

    private:
      const vector<Primitive>& edge_slopes_;
    } minmoder(edge_slopes);

    return serial_generate(minmoder);
  }
}

VanLeer::VanLeer(void) {}

vector<std::pair<Primitive,Primitive> > 
VanLeer::interpolateAll
(const HydroSnapshot& hs,
 double /*dt*/) const
{
  const class Interpolator: 
  public Index2Member<std::pair<Primitive,Primitive> >
  {
  public:

    Interpolator(const HydroSnapshot& hs_i):
      hs_(hs_i),
      cell_slopes_(calc_cell_slopes(hs_i)) {}

    size_t getLength(void) const
    {
      return hs_.edges.size()-2;
    }

    std::pair<Primitive, Primitive> operator()(size_t i) const
    {
      return std::pair<Primitive, Primitive>
	(hs_.cells[i]+
	 0.5*(hs_.edges[i+1]-hs_.edges[i])*cell_slopes_[i],
	 hs_.cells[i+1]-
	 0.5*(hs_.edges[i+2]-hs_.edges[i+1])*cell_slopes_[i+1]);
    }

  private:
    const HydroSnapshot& hs_;
    const vector<Primitive> cell_slopes_;
  } interpolator(hs);

  return serial_generate(interpolator);
}
