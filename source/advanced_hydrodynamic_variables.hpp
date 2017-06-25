/*! \file advanced_hydrodynamic_variables.hpp
  \author Almog Yalinewich
  \brief Functions that require both hydrodyanmic_variables and equation_of_state
  \details One would naively think that these functions should be included 
  in hydrodynamic_variables, but I decided to use another file to avoid cyclic dependence
*/

#ifndef ADVANCED_HYDRODYNAMIC_VARIABLES
#define ADVANCED_HYDRODYNAMIC_VARIABLES 1

#include "hydrodynamic_variables.hpp"
#include "equation_of_state.hpp"

/*! \brief Converts primitives to conserved variables
  \param p Primitive variables
  \param eos Equation of state
  \return Conserved variables
 */
Conserved Primitive2Conserved(Primitive const& p, EquationOfState const& eos);

/*! \brief Calculates new conserved variables from the primitive variables
  \param p Primitive variables
  \param eos Equation of state
  \return new conserved variables
 */
NewConserved primitive_to_new_conserved(const Primitive& p, const EquationOfState& eos);

/*! \brief Converts old conserved variables to new
  \param c Old conserved variables
  \return New conserved variables
 */
NewConserved old_to_new_conserved(const Conserved& c);

/*! \brief Converts new conserved variables to old conserved variables
  \param c New conserved variables
  \return Old conserved variables
 */
Conserved new_to_old_conserved(const NewConserved& c);

/*! \brief Converts primitives to flux
  \param hs Primitive variables
  \param eos Equation of state
  \return Flux
 */
Conserved Primitive2Flux(Primitive const& hs, EquationOfState const& eos);

/*! \brief Converts primitives to flux
  \param hs Primitive variables
  \param eos Equation of state
  \return Flux
 */
NewConserved primitive_to_flux(Primitive const& hs, EquationOfState const& eos);

/*! \brief Converts primitives to conserved per volume (instead of per particle)
  \param p Primitive variables
  \param eos Equation of state
  \return Conserved variables (per volume)
 */
Conserved Primitive2Conserved_pv(Primitive const& p,
				 const EquationOfState& eos);

/*! \brief Converts primitives to conserved per volume (instead of per particle)
  \param p Primitive variables
  \param eos Equation of state
  \return Conserved variables (per volume)
 */
NewConserved primitive_to_conserved_pv(Primitive const& p,
				       const EquationOfState& eos);

/*! \brief conserts primitives to fluxes (per volume)
  \param p Primitive variables
  \param eos Equation of state
  \return Fluxes per volume
 */
Conserved Primitive2Flux_pv(Primitive const& p,
			    const EquationOfState& eos);

/*! \brief conserts primitives to fluxes (per volume)
  \param p Primitive variables
  \param eos Equation of state
  \return Fluxes per volume
 */
NewConserved primitive_to_flux_pv(const Primitive& p,
				  const EquationOfState& eos);

#endif // ADVANCED_HYDRODYNAMIC_VARIABLES
