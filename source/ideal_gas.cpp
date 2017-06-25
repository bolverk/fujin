#include <limits>
#include <cmath>
#include <cassert>
#include "utilities.hpp"
#include "ideal_gas.hpp"
#include "enthalpy_equation.hpp"
#include "advanced_hydrodynamic_variables.hpp"
#include "universal_error.hpp"
#include "aimm_recovery.hpp"
#include <assert.h>

// Checks that the energy momentum and momentum densities are physically feasible

namespace {
  /*! \brief Checks if the input variables are physically valid
    \param positive Sum of energy and momentum
    \param negative Difference of energy and momentum
    \return True if em and sm are physical
   */
  bool are_valid_parameters(double positive, double negative)
  {
    const bool cond1 = positive+negative>2;
    const bool cond2 = positive*negative>1;
    return (cond1&&cond2);
  }
}

// Enthalpy equation

// Solves the enthalpy equation in the limit of the cold - fast approximation
namespace {

  /*! \brief Calculates the enthalpy in the case of cold, fast fluid
    \param positive Sum of energy and momentum
    \param negative Difference of energy and momentum
    \param g Adiabatic index
    \return Enthalpy
   */
  double calc_dh_cold_fast(double positive, double negative, double g)
  {
    const double em = 0.5*(positive+negative);
    const double sm = 0.5*(positive-negative);
    const double res = g*sqrt(pow(sm,2)+1)*
      (positive*negative-1)/
      (1+g-em*g/sqrt(1+pow(sm,2)))/
      (em+sqrt(pow(sm,2)+1));
    return res;
  }

  /*! \brief Determines whether to use the cold fast approximation instead of the exact solution
    \param positive Energy plus momentum
    \param negative Energy minum momentum
    \param g Adiabatic index
    \return True if cold fast approximation should be used
   */
  bool cold_fast_fluid_criterion(double positive, double negative, double g)
  {
    const double em = 0.5*(positive+negative);
    const double sm = 0.5*(positive-negative);
    return (em>100) && fabs(em-fabs(sm))<1 && 
      (fabs(calc_dh_cold_fast(positive,negative,g))<0.05);
  }
}

double conserved2enthalpy_diff(double positive, double negative, double g)
{
  // Main
  const double em = 0.5*(positive + negative);
  const double sm = 0.5*(positive - negative);
  if(cold_fast_fluid_criterion(positive, negative, g))
    return calc_dh_cold_fast(positive, negative, g);
  else if(fabs(sm/em)<1e-6)
    return g*(em-1);
  else{

    // Initial guess based on approximation for large enthalpy
    const double dhguess = g*(em-1);
    // Enthalpy equation
    const EnthalpyEquation ee(em, sm, g);

    const SVRelStep sc(1e-12);
    const double res = NRSafe(ee, 0,dhguess*1.01,sc);

    assert(res>=0);

    return res;
  }
}

// Public

IdealGas::IdealGas(double ig):
  g(ig) {}

double IdealGas::dp2e(double d, double p) const
{
  return p/(g-1)+d;
}

double IdealGas::dp2ba(double d, double p) const
{
  /*
  if(std::numeric_limits<double>::epsilon()>d && 
     std::numeric_limits<double>::epsilon()>p)
    return 0;
  */
  const double e = dp2e(d, p);
  return sqrt(g*p/(e+p));
}

Primitive IdealGas::Conserved2Primitive
(Primitive const& old, NewConserved const& c) const
{
  assert(c.mass>0);
  if(are_valid_parameters(c.positive,c.negative)){
    /*
    const double dh2 = conserved2enthalpy_diff
      (c.positive, c.negative, g);
    const double h2 = 1 + dh2;
    const double w2 = 0.5*(c.positive-c.negative)/h2;
    const double lf2 = celerity2lorentz_factor(w2);
    const double d2 = c.mass/lf2;
    const double p2 = dh2*(g-1)*d2/g;
    assert(p2>0);
    const Primitive ref(d2,p2,w2);
    */

    const double p = calc_pressure(c,g);
    const double p2D = p/c.mass;
    const double w = 0.5*(c.positive-c.negative)/
      sqrt((c.positive+p2D)*(c.negative+p2D));
    const double d = c.mass/celerity2lorentz_factor(w);

    /*
    assert(fabs(d-d2)/(d+d2)<1e-6 &&
	   fabs(p-p2)/(p+p2)<1e-6 &&
	   fabs(w-w2)<1e-6 && 
	   "new method does not reproduce old results");
    */

    return Primitive(d,p,w);
  }
  else{
    const double s = old.Pressure/pow(old.Density,g);
    const double w = 0.5*(c.positive - c.negative);
    const double lf = celerity2lorentz_factor(w);
    const double d = c.mass/lf;
    const double p = s*pow(d,g);
    return Primitive(d,p,w);
  }
 
  /* 
  const Conserved temp = Primitive2Conserved(Primitive(d,p,w),*this);
  const double aux1 = fabs(temp.Energy-c.Energy)/(temp.Energy+c.Energy);
  assert(aux1<1e-8 && "Failed to reproduce conserved in ideal gas");
  */
}

double IdealGas::getAdiabaticIndex(void) const
{
  return g;
}
