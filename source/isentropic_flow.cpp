#include <cmath>
#include "isentropic_flow.hpp"
#include "ideal_gas.hpp"

ConstEntropy::ConstEntropy
(double k, double g,
 SpatialDistribution const& density):
  k_(k), g_(g), density_(density) {}

double ConstEntropy::operator()(double x) const
{
  const double d = density_(x);
  return k_*pow(d,g_);
}

// Constant riemann invatiant

double calc_riemann_invariant
(double g, double d, double p,
 double v, double s)
{
  const double cs = IdealGas(g).dp2ba(d,p);
  return atanh(v)+s*(2/sqrt(g-1))*
    atanh(cs/sqrt(g-1));
}

ConstRiemannInv::ConstRiemannInv
(double jm, double g,
 SpatialDistribution const& density,
 SpatialDistribution const& pressure):
  jm_(jm), g_(g),
  density_(density),
  pressure_(pressure) {}

double ConstRiemannInv::operator()(double x) const
{
  const double d = density_(x);
  const double p = pressure_(x);
  const double aux = jm_ - 
    calc_riemann_invariant(g_,d,p,0,-1);
  return tanh(aux);
}
