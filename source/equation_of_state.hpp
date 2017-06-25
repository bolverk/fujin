/*! \file equation_of_state.hpp
  \brief Equation of state module
  \author Almog Yalinewich
*/

#ifndef EQUATION_OF_STATE_HPP
#define EQUATION_OF_STATE_HPP 1

#include "hydrodynamic_variables.hpp"

//! \brief Base class for equations of state
class EquationOfState
{
public:
  /*! \brief Calculates the energy density
    \param d Density
    \param p Pressure
    \return Energy per unit volume
   */
  virtual double dp2e(double d, double p) const = 0;
  /*! \brief Calculates the dimensionless speed of sound
    \param d Density
    \param p Pressure
    \return Dimensionless speed of sound
   */
  virtual double dp2ba(double d, double p) const = 0;  

  /*! \brief Converts conserved variables to primitives
    \param old Value of the primitives in the previous cycle
    \param c Conserved variables
    \return Primitive variables
   */
  virtual Primitive Conserved2Primitive
  (Primitive const& old, NewConserved const& c) const = 0;

  virtual ~EquationOfState(void);
};

#endif
