/*! \file linear_rs.hpp
  \brief Linear Riemann solver
  \author Almog Yalinewich
 */

#ifndef LINEAR_RS_HPP
#define LINEAR_RS_HPP 1

#include "riemann_solver.hpp"
#include "equation_of_state.hpp"

//! \brief Linear Riemann solver
class LinearRS: public RiemannSolver
{
public:

  /*! \brief Class constructor
    \param eos Equation of state
   */
  explicit LinearRS(const EquationOfState& eos);

  RiemannSolution operator()
  (const Primitive& left, const Primitive& right) const override;
  
private:

  //! \brief Equation of state
  const EquationOfState& eos_;
};

#endif // LINEAR_RS_HPP
