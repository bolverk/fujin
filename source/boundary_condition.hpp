/*! \brief Boundary condition
  \file boundary_condition.hpp
  \author Almog Yalinewich
 */

#ifndef BOUNDARY_CONDITION_HPP
#define BOUNDARY_CONDITION_HPP 1

#include <vector>
#include "riemann_solver.hpp"

//! \brief Base class for boundary condition
template<template<class> class CP>
class BoundaryCondition
{
public:
  /*! \brief Calculates the Riemann solution at the edges
    \param side False for left, true for right
    \param cells Primitive variables
    \return Riemann solution (pressure and velocity)
   */
  virtual RiemannSolution
  operator()(bool side, const CP<Primitive>& cells) const = 0;

  virtual ~BoundaryCondition(void) {}
};

#endif // BOUNDARY_CONDITION_HPP
