#ifndef VAN_LEER_HPP
#define VAN_LEER_HPP 1

#include "spatial_reconstruction.hpp"
#include "utilities.hpp"

namespace van_leer
{

  /*! \brief returns the sign of a number
    \param x Number
    \return Sign
  */
  static double sgn(double x)
  {
    return (x>0) - (x<0);
  }

  /*! \brief Minmod
    \param x First argument
    \param y Second argument
    \return Minmod of x and y
  */
  static inline double minmod(double x, double y)
  {
    return 0.5*(sgn(x)+sgn(y))*
      min(fabs(x),fabs(y));
  }

  /*! \brief Applies minmod filter to all members
    \param p1 First set of primitive variables
    \param p2 Second set of primitive variables
    \return Minmod filtered primitve variables 
   */
  static inline Primitive minmod(Primitive const& p1,
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
  template<template<class> class CE, template<class> class CP>
  vector<Primitive> calc_edge_slopes
  (const NewHydroSnapshot<CE, CP>& hs)
  {
    const class DerivativeInterface: public Index2Member<Primitive>
    {
    public:

      explicit DerivativeInterface
      (const NewHydroSnapshot<CE, CP>& hs_i):
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
      const NewHydroSnapshot<CE, CP>& hs_;
    } deriv_interface(hs);
    return serial_generate(deriv_interface);
  }

  /*! \brief Calculates the slopes on the computational cells
    \param hs Hydrodynamic snapshot
    \return List of slopes
   */
  template<template<class> class CE, template<class> class CP>
  vector<Primitive> calc_cell_slopes
  (const NewHydroSnapshot<CP, CE>& hs)
  {
    const vector<Primitive> edge_slopes = calc_edge_slopes(hs);
    const class Minmoder: public Index2Member<Primitive>
    {
    public:

     explicit Minmoder(const vector<Primitive>& edge_slopes_i):
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

/*! \brief Spatial reconstruction based on Van Leer's method
  \details The description of the method is available in <a href="http://www.sciencedirect.com/science/article/pii/S0021999197957041">B. Van Leer, 'Towards the Ultimate Conservative Difference Scheme', J. Comp. Phys. 135, 229-248 (1977) </a>
*/
template<template<class> class CE, template<class> class CP>
class VanLeer: public SpatialReconstruction<CE, CP>
{
public:

  /*! \brief Class constructor
   */
  VanLeer(void) {}

  vector<std::pair<Primitive,Primitive> > interpolateAll
  (const NewHydroSnapshot<CE, CP>& hs,
   double /*dt*/) const
  {
    const class Interpolator: 
    public Index2Member<std::pair<Primitive,Primitive> >
    {
    public:

      explicit Interpolator
      (const NewHydroSnapshot<CE, CP>& hs_i):
	hs_(hs_i),
	cell_slopes_(van_leer::calc_cell_slopes(hs_i)) {}

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
      const NewHydroSnapshot<CE, CP>& hs_;
      const vector<Primitive> cell_slopes_;
    } interpolator(hs);

    return serial_generate(interpolator);
  }
};

#endif // VAN_LEER_HPP
