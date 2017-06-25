#ifndef FREE_FLOW_HPP
#define FREE_FLOW_HPP 1

#include "boundary_condition.hpp"
#include "riemann_solver.hpp"

/*! \brief Free flow boundary (continuous) conditions
 */
class FreeFlow: public BoundaryCondition
{
  
public:

  RiemannSolution CalcRS
  (size_t idx, vector<Primitive> const& cells) const;
};

#endif // FREE_FLOW_HPP
