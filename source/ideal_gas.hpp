/*! \file ideal_gas.hpp
  \brief Ideal gas equation of state
  \author Almog Yalinewich
 */

#ifndef IDEAL_GAS_HPP
#define IDEAL_GAS_HPP 1

#include "equation_of_state.hpp"
#include "hydrodynamic_variables.hpp"

//! \brief Ideal gas equation of state
class IdealGas: public EquationOfState
{
private:
  //! \brief Adiabatic index
  const double g;
public:
  /*! \brief Class constructor
    \param ig Adiabatic index
   */
  explicit IdealGas(double ig);

  double dp2e(double d, double p) const;

  double dp2ba(double d, double p) const;

  Primitive Conserved2Primitive
  (Primitive const& old, NewConserved const& c) const;

  /*! \brief Returns the value of the adiabatic index
   */
  double getAdiabaticIndex(void) const;
};

/*! \brief Calculates the enthalpy
  \param em Energy
  \param sm Momentum
  \param g Adiabatic index
  \return Enthalpy
 */
double conserved2enthalpy_diff(double em, double sm, double g);

#endif // IDEAL_GAS_HPP
