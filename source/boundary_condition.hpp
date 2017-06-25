/*! \brief Boundary condition
  \file boundary_condition.hpp
  \author Almog Yalinewich
 */

#ifndef BOUNDARY_CONDITION_HPP
#define BOUNDARY_CONDITION_HPP 1

#include <vector>
#include "riemann_solver.hpp"

//! \brief Base class for boundary condition
class BoundaryCondition
{
public:
  /*! \brief Calculates the Riemann solution at the edges
    \param idx Vertex index
    \param cells Primitive variables
    \return Riemann solution (pressure and velocity)
   */
  virtual RiemannSolution
  CalcRS(size_t idx, vector<Primitive> const& cells) const = 0;

  virtual ~BoundaryCondition(void);

  //! \brief Invalid index error
  class InvalidIndexError
  {
  public:
    //! \brief Invalid index
    unsigned int Index;
  };
};

#endif // BOUNDARY_CONDITION_HPP
