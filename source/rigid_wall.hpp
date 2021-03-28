#ifndef RIGID_WALL_HPP
#define RIGID_WALL_HPP 1

#include "boundary_condition.hpp"
#include "riemann_solver.hpp"

/*! \brief Rigid wall boundary conditions
 */
template<template<class> class CP>
class RigidWall: public BoundaryCondition<CP>
{

private:

  Primitive invert_celerity(const Primitive& p) const
  {
    Primitive res = p;
    res.Celerity = -res.Celerity;
    return res;
  }
  
public:

  /*! \brief Class constructor
    \param rs Pointer to Riemann solver
   */
  explicit RigidWall(RiemannSolver const& rs):
    rs_(rs) {}

  RiemannSolution operator()
  (bool side, CP<Primitive> const& cells) const override
  {
    return side ? rs_(cells.back(), invert_celerity(cells.back())) :
      rs_(invert_celerity(cells.front()), cells.front());
  }

private:

  //! \brief Pointer to Riemann solver
  RiemannSolver const& rs_;
};

#endif // RIGID_WALL_HPP
