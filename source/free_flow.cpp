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

RiemannSolution FreeFlow::operator()
(bool side, vector<Primitive> const& cells) const
{
  return side ? Primitive2RiemannSolution(cells.back()) :
    Primitive2RiemannSolution(cells.front());
}
