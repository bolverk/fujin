/*! \file aimm_recovery.hpp
  \author Almog Yalinewich
  \brief Calculates the pressure from the conserved variables
  \details Based on appendix C in http://adsabs.harvard.edu/abs/1999ApJS..122..151A
 */

#ifndef AIMM_RECOVRY_HPP
#define AIMM_RECOVRY_HPP 1

#include "hydrodynamic_variables.hpp"

/*! \brief Calculates the pressure
  \param c Conserved variables
  \param g Adiabatic index
  \return Pressure
 */
double calc_pressure(const NewConserved& c,
		     double g);

#endif 
