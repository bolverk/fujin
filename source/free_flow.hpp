#ifndef FREE_FLOW_HPP
#define FREE_FLOW_HPP 1

#include "boundary_condition.hpp"
#include "riemann_solver.hpp"

/*! \brief Free flow boundary (continuous) conditions
 */
template<template<class> class CP>
class FreeFlow: public BoundaryCondition<CP>
{
private:
  
  RiemannSolution Primitive2RiemannSolution(const Primitive& p) const
  {
    return RiemannSolution(p.Pressure,p.Celerity);
  }
  
public:

  RiemannSolution operator()
  (bool side,
   const CP<Primitive>& cells) const override
  {
    return side ? Primitive2RiemannSolution(cells.back()) :
      Primitive2RiemannSolution(cells.front());
  }
};

#endif // FREE_FLOW_HPP
