//! \file imgrs.hpp
//! \brief Exact Riemann solver for relativistic ideal gas
//! \author Almog Yalinewich

#ifndef IMGRS_HPP
#define IMGRS_HPP 1

#include "riemann_solver.hpp"
#include "utilities.hpp"

//! \brief Exact Riemann solver for an ideal gas
class IdealGasRiemannSolver : public RiemannSolver
{
private:
  //! \brief Adiabatic index
  double g;
public:
  /*! \brief Class constructor
    \param ig Adiabatic index
  */
  explicit IdealGasRiemannSolver(double ig);
  /*! \brief Solve the Riemann problem
    \param left Hydrodynamic state on the left side
    \param right Hydrodynamic state on the right side
    \return Pressure and velocity at the interface
   */
  RiemannSolution operator()(Primitive const& left, 
			     Primitive const& right) const;
};

#endif
