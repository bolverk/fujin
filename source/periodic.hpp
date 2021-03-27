#ifndef PERIODIC_HPP
#define PERIODIC_HPP 1

#include "boundary_condition.hpp"
#include "riemann_solver.hpp"

//! \brief Periodic boundary conditions
template<template<class> class CP>
class Periodic: public BoundaryCondition<CP>
{
public:

  /*! \brief Class constructor
    \param rs Riemann solver
   */
  explicit Periodic(RiemannSolver const& rs):
    rs_(rs) {}

  RiemannSolution operator()
  (bool side, const CP<Primitive>& cells) const override
  {
    return side ? rs_(cells.front(),cells.back()) :
     rs_(cells.back(),cells.front());
  }
  
private:

  //! \brief Riemann solver  
  const RiemannSolver& rs_;
};

# endif // PERIODIC_HPP
