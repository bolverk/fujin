#ifndef ISENTROPIC_FLOW_HPP
#define ISENTROPIC_FLOW_HPP 1

#if 0
#include "spatial_distribution.hpp"

//! \brief Pressure distribution that keeps the entropy constant
class ConstEntropy: public SpatialDistribution
{
public:

  /*! \brief Class constructor
    \param k Prefactor
    \param g Adiabatic index
    \param density Density distribution
   */
  ConstEntropy(double k, double g,
	       SpatialDistribution const& density);

  double operator()(double x) const;

private:

  //! \brief Prefactor
  const double k_;

  //! \brief Adiabatic index
  const double g_;

  //! \brief Density distribution
  SpatialDistribution const& density_;
};

/*! \brief Calculates the Riemann invariant
  \param g Adiabatic index
  \param d Density
  \param p Pressure
  \param v Velocity
  \param s direction
  \return Riemann invariant
 */
double calc_riemann_invariant
(double g, double d, double p,
 double v, double s);

//! \brief Velocity distribution that keeps the riemann invariant constant
class ConstRiemannInv: public SpatialDistribution
{
public:

  /*! \brief Class constructor
    \param jm Riemann invariant
    \param g Adiabatic index
    \param density Density distribution
    \param pressure Pressure distribution
   */
  ConstRiemannInv
  (double jm, double g,
   SpatialDistribution const& density,
   SpatialDistribution const& pressure);

  double operator()(double x) const;

private:

  //! \brief Riemann invariant
  const double jm_;

  //! \brief Adiabatic index
  const double g_;

  //! \brief Density distribution
  SpatialDistribution const& density_;

  //! \brief Pressure distribution
  SpatialDistribution const& pressure_;
};
#endif // 0
#endif // ISENTROPIC_FLOW_HPP
