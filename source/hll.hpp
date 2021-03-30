/*! \file hll.hpp
  \brief HLL Riemann solver
  \details Based on \n 
  <a href="http://adsabs.harvard.edu/full/2005MNRAS.364..126M">A. Mignone & G. Bodo, "An HLLC Riemann solver for relativistic flows - I. Hydrodynamics", Roy. Astr. Soc. vol 364, issue 1, pp 126-136</a> 
  \author Almog Yalinewich
 */

#ifndef HLL_HPP
#define HLL_HPP 1

#include "riemann_solver.hpp"
#include "equation_of_state.hpp"

//! \brief HLL Riemann solver
class HLL: public RiemannSolver
{
private:

public:
  /*! \brief Class constructor
    \param eos Equation of state
   */
  explicit HLL(const EquationOfState& eos);

  /*! \brief Solve the Riemann problem
    \param left Hydrodynamic state on the left side
    \param right Hydrodynamic state on the right side
    \return Pressure and velocity at the interface
   */
  RiemannSolution operator()(Primitive const& left, 
			     Primitive const& right) const override;

private:

  //! \brief Equation of state  
  const EquationOfState& eos_;
};

#endif
