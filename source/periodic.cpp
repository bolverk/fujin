#include <cassert>
#include "periodic.hpp"

Periodic::Periodic(RiemannSolver const& rs):
  rs_(rs) {}

RiemannSolution Periodic::CalcRS(size_t idx,
				 vector<Primitive> const& cells) const
{
  assert(idx==0 || cells.size()==idx);
  return idx==0 ? rs_(cells.back(),cells.front()) :
    rs_(cells.front(),cells.back());
}
