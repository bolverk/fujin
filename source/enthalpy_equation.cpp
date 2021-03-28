#include "enthalpy_equation.hpp"

namespace {
  /*! \brief Calcuates the coefficients
    \param em Energy per particle (lab frame)
    \param sm Momentum per particle (lab frame)
    \param g Adiabatic index
    \return List of coefficients
   */
  /*
  vector<double> calc_coefs(double em, double sm, double g)
  {
    vector<double> res;
    res.push_back(1);
    res.push_back(2*(g+1));
    const double de2 = (em-sm)*(em+sm);
    const double aux = pow(sm,2)+1;
    res.push_back(1+2*(1+aux)*g-pow(g,2)*(-2+aux+de2));
    res.push_back(2*g*(aux+g-g*de2));
    res.push_back(-aux*pow(g,2)*(de2-1));
    return res;
  }
  */

  quartic<double> calc_coefs(double em, double sm, double g)
  {
    quartic<double> res;
    res[0] = 1;
    res[1] = 2*(g+1);
    const double de2 = (em-sm)*(em+sm);
    const double aux = pow(sm,2)+1;
    res[2] = 1+2*(1+aux)*g-pow(g,2)*(-2+aux+de2);
    res[3] = 2*g*(aux+g-g*de2);
    res[4] = -aux*pow(g,2)*(de2-1);
    return res;
  }
}

EnthalpyEquation::EnthalpyEquation(double em, double sm, double g):
  poly_(calc_coefs(em,sm,g)) {}

double EnthalpyEquation::operator()(double dh) const
{
  return poly_(dh);
}

double EnthalpyEquation::diff(double dh) const
{
  return poly_.diff(dh);
}
