#include <cassert>
#include "rigid_wall.hpp"
#include "advanced_hydrodynamic_variables.hpp"

RigidWall::RigidWall(RiemannSolver const& rs):
  rs_(rs) {}

namespace {

  /*! \brief Reverses the celerity
    \param p Primitive variable
    \return Primitive variable with reversed celerity
   */
  Primitive invert_celerity(Primitive const& p)
  {
    Primitive res = p;
    res.Celerity = -res.Celerity;
    return res;
  }
}

RiemannSolution RigidWall::CalcRS
(size_t idx, vector<Primitive> const& cells) const
{
  assert(idx==0 || idx==cells.size());
  return idx==0 ? rs_(invert_celerity(cells.front()), cells.front()) :
    rs_(cells.back(), invert_celerity(cells.back()));
}
