#ifndef COLLELA_HPP
#define COLLELA_HPP 1

#include "spatial_distribution.hpp"

//! \brief Spatial distribution attributed to P. Collela
class Collela: public SpatialDistribution
{
public:

  /*! \brief Class constructor
    \param ref Vertical offset
    \param alpha Amplitude
    \param l Width
    \param x0 Horizontal offset
   */
  Collela(double ref, double alpha, 
	  double l, double x0);

  double operator()(double x) const;

private:

  //! \brief vertical offset
  const double ref_;

  //! \brief Amplitude
  const double a_;

  //! \brief Width
  const double l_;

  //! \brief Horizontal offset
  const double x0_;
};

#endif // COLLELA_HPP
