#ifndef RIGID_WALL_HPP
#define RIGID_WALL_HPP 1

#include "boundary_condition.hpp"
#include "riemann_solver.hpp"

/*! \brief Rigid wall boundary conditions
 */
class RigidWall: public BoundaryCondition
{
  
public:

  /*! \brief Class constructor
    \param rs Pointer to Riemann solver
   */
  RigidWall(RiemannSolver const& rs);

  RiemannSolution CalcRS
  (size_t idx, vector<Primitive> const& cells) const;

private:

  //! \brief Pointer to Riemann solver
  RiemannSolver const& rs_;
};

#endif // RIGID_WALL_HPP
