#ifndef PERIODIC_HPP
#define PERIODIC_HPP 1

#include "boundary_condition.hpp"
#include "riemann_solver.hpp"

//! \brief Periodic boundary conditions
class Periodic: public BoundaryCondition
{
public:

  /*! \brief Class constructor
    \param rs Riemann solver
   */
  explicit Periodic(RiemannSolver const& rs);

  RiemannSolution CalcRS
  (size_t idx, vector<Primitive> const& cells) const;

private:

  //! \brief Riemann solver  
  RiemannSolver const& rs_;
};

# endif // PERIODIC_HPP
