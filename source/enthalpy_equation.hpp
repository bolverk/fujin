/*! \file enthalpy_equation.hpp
  \brief Enthalpy equation for ideal gas
  \author Almog Yalinewich
*/

#ifndef ENTHALPY_EQUATION_HPP
#define ENTHALPY_EQUATION_HPP 1

#include "utilities.hpp"
#include "polynomial.hpp"

template<class T> using quartic = array<double, 5>;

//! \brief Enthalpy equation
class EnthalpyEquation: public SVDifferentiable
{
public:

  /*! \brief Class constructor
    \param em Energy density
    \param sm Momentum density
    \param g Adiabatic index
   */
  EnthalpyEquation(double em, double sm, double g);

  double operator()(double h) const;

  double diff(double h) const;

private:

  //! \brief Left hand side of enthalpy equation
  const Polynomial<quartic> poly_;
};

#endif // ENTHALPY_EQUATION_HPP
