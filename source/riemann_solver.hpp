/*! \file riemann_solver.hpp
  \author Almog Yalinewich
  \brief Abstract class for a Riemann solver
 */

#ifndef RIEMANN_SOLVER_HPP
#define RIEMANN_SOLVER_HPP 1

#include "hydrodynamic_variables.hpp"

//! \brief Solution to the Riemann problem
struct RiemannSolution
{
  /*! \brief Void constructor
   */
  RiemannSolution(void);

  /*! \brief Class constructor
    \param p Pressure
    \param w Celerity
   */
  RiemannSolution(double p, double w);

  //! \brief Pressure
  double Pressure;

  //! \brief Dimensionless velocity
  double Celerity;
};

//! \brief Base class for a Riemann solver
class RiemannSolver
{
public:
  /*! \brief Solve the Riemann problem
    \param left Hydrodynamic state on the left side
    \param right Hydrodynamic state on the right side
    \return Pressure and velocity at the interface
  */
  virtual RiemannSolution operator()
  (Primitive const& left, 
   Primitive const& right) const = 0;

  virtual ~RiemannSolver(void);
};

#endif
