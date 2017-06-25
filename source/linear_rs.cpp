#include <cmath>
#include "linear_rs.hpp"

LinearRS::LinearRS(const EquationOfState& eos):
  eos_(eos) {}

namespace {

  /*! \brief Calculates the weights for the arithmetic average
    \param cell Primitive variables
    \param eos Equation of state
    \return Weight
   */
  double calc_weight(const Primitive& cell,
		     const EquationOfState& eos)
  {
    const double h = cell.Pressure + eos.dp2e(cell.Density,
					      cell.Pressure);
    const double ba = eos.dp2ba(cell.Density, cell.Pressure);
    const double w = cell.Celerity;
    return h*ba/sqrt(1+pow(w,2));
  }
}

RiemannSolution LinearRS::operator()(const Primitive& left,
				     const Primitive& right) const
{
  const double weight_l = calc_weight(left, eos_);
  const double weight_r = calc_weight(right, eos_);
  const double w_s = 
    (left.Celerity*weight_l+right.Celerity*weight_r)/(weight_r+weight_l)+
    (left.Pressure - right.Pressure)/(weight_l+weight_r);
  const double iweight_l = 1./weight_l;
  const double iweight_r = 1./weight_r;
  const double p_s =
    (left.Pressure*iweight_l + right.Pressure*iweight_r)/(iweight_r+iweight_l)+
    (left.Celerity - right.Celerity)/(iweight_l+iweight_r);
  return RiemannSolution(p_s, w_s);
}
