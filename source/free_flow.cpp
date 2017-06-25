#include "free_flow.hpp"
#include <cassert>

namespace {

  /*! \brief Creates an instance of a Riemann solution
    \param p Primitive
    \return Riemann solution
   */
  RiemannSolution Primitive2RiemannSolution(Primitive const& p)
  {
    return RiemannSolution(p.Pressure,p.Celerity);
  }
}

RiemannSolution FreeFlow::CalcRS
(size_t idx, vector<Primitive> const& cells) const
{
  assert(idx==0 || idx==cells.size());
  return idx==0 ? Primitive2RiemannSolution(cells.front()) :
    Primitive2RiemannSolution(cells.back());
}
