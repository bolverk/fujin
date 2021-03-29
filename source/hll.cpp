/*! \file hll.cpp
  \brief Implementation of hll.hpp
  \author Almog Yalinewich
*/

#include <cmath>
#include "hll.hpp"
#include "utilities.hpp"
#include "advanced_hydrodynamic_variables.hpp"
#include "linear_rs.hpp"

using std::min;
using std::max;
using std::abs;

HLL::HLL(const EquationOfState& eos):
  eos_(eos) {}

namespace {
  /*! \brief Calculates the right wave velocity
    \param left Flow variables on the left side
    \param right Flow variables on the right side
    \param eos Equation of state
    \return Left wave velocity
  */
  /*
  double lambda_l(Primitive const& left, 
		  Primitive const& right,
		  const EquationOfState& eos)
  {
    const double left_velocity = celerity2velocity(left.Celerity);
    const double right_velocity = celerity2velocity(right.Celerity);
    const double left_sound_speed = eos.dp2ba(left.Density,
					       left.Pressure);
    const double right_sound_speed = eos.dp2ba(right.Density,
						right.Pressure);
    return min(RelVelAdd(left_velocity,-left_sound_speed),
	       RelVelAdd(right_velocity,-right_sound_speed));
  }
  */

  /*! \brief Calculates the right wave velocity
    \param left Flow variables on the left side
    \param right Flow variables on the right side
    \param eos Equation of state
    \return Right wave velocity
  */
  /*
  double lambda_r(Primitive const& left,
		  Primitive const& right,
		  const EquationOfState& eos)
  {
    const double left_velocity = celerity2velocity(left.Celerity);
    const double right_velocity = celerity2velocity(right.Celerity);
    const double left_sound_speed = eos.dp2ba(left.Density,
					      left.Pressure);
    const double right_sound_Speed = eos.dp2ba(right.Density,
					       right.Pressure);
    return max(RelVelAdd(left_velocity,left_sound_speed),
	       RelVelAdd(right_velocity,right_sound_Speed));
  }
  */

  /*! \brief Calculates the HLL conserved variables
    \param left Flow variables on the left side
    \param right Flow variables on the right side
    \param eos Equation of state
    \return HLL conserved variables
  */
  /*
  NewConserved CalcUhll(Primitive const& left, 
		     Primitive const& right,
		     const EquationOfState& eos)
  {
    const double l_l = lambda_l(left, right, eos);
    const double l_r = lambda_r(left, right, eos);
    const NewConserved u_l = 
      primitive_to_conserved_pv(left,eos);
    const NewConserved u_r = 
      primitive_to_conserved_pv(right,eos);
    const NewConserved f_l = 
      primitive_to_flux_pv(left,eos);
    const NewConserved f_r =
      primitive_to_flux_pv(right,eos);

    return (1./(l_r-l_l))*(l_r*u_r-l_l*u_l+f_l-f_r);
  }
  */

  /*! \brief Calculates the HLL flux
    \param left Flow variables on the left side
    \param right Flow variables on the right side
    \param eos Equation of state
    \return HLL fluxes
  */
  /*
  NewConserved CalcFhll(Primitive const& left, 
			Primitive const& right,
			const EquationOfState& eos)
  {
    const double l_l = lambda_l(left, right, eos);
    const double l_r = lambda_r(left, right, eos);
    const NewConserved u_l = 
      primitive_to_conserved_pv(left,eos);
    const NewConserved u_r = 
      primitive_to_conserved_pv(right,eos);
    const NewConserved f_l = 
      primitive_to_flux_pv(left,eos);
    const NewConserved f_r =
      primitive_to_flux_pv(right,eos);
    return (1./(l_r-l_l))*
      (l_r*f_l-l_l*f_r+l_r*l_l*(u_r-u_l));
  }
  */

  std::pair<double,double> calc_lambdas(const Primitive& hv,
					const EquationOfState& eos)
  {
    const double ba = eos.dp2ba(hv.Density, hv.Pressure);
    const double lf = celerity2lorentz_factor(hv.Celerity);
    const double sigma = pow(ba/lf,2)/(1-pow(ba,2));
    const double v = hv.Celerity/lf;
    const double discriminant = 
      sqrt(sigma*(1-pow(v,2)+sigma));
    return std::pair<double,double>
      ((v+discriminant)/(1+sigma),
       (v-discriminant)/(1+sigma));
  }

  std::pair<double, double> estimate_wave_velocities
  (const Primitive& left, 
   const Primitive& right,
   const EquationOfState& eos)
  {
    const std::pair<double, double> left_velocities =
      calc_lambdas(left, eos);
    const std::pair<double, double> right_velocities =
      calc_lambdas(right, eos);
    return std::pair<double, double>
      (min(left_velocities.second,
	   right_velocities.second),
       max(left_velocities.first,
	   right_velocities.first));
  }

  std::pair<NewConserved,NewConserved> calc_uhll_fhll
  (const Primitive& left,
   const Primitive& right,
   const EquationOfState& eos)
  {
    const std::pair<double, double> wave_speeds =
      estimate_wave_velocities(left, right, eos);
    const double l_l = wave_speeds.first;
    const double l_r = wave_speeds.second;
    /*
    const double l_l = lambda_l(left, right, eos);
    const double l_r = lambda_r(left, right, eos);
    */
    const NewConserved u_l = 
      primitive_to_conserved_pv(left,eos);
    const NewConserved u_r = 
      primitive_to_conserved_pv(right,eos);
    const NewConserved f_l = 
      primitive_to_flux_pv(left,eos);
    const NewConserved f_r =
      primitive_to_flux_pv(right,eos);
    return std::pair<NewConserved,NewConserved>
      ((1./(l_r-l_l))*(l_r*u_r-l_l*u_l+f_l-f_r),
       (1./(l_r-l_l))*
       (l_r*f_l-l_l*f_r+l_r*l_l*(u_r-u_l)));
  }

  /*! \brief Calculates the celerity at the interface
    \param uhll HLL conserved variables
    \param fhll HLL fluxes
    \return Celerity
   */
  double calc_hll_ws(const NewConserved& uhll,
		     const NewConserved& fhll)
  {
    const double thres = 1e-6;
    /*
    const double um = uhll.Momentum;
    const double ue = uhll.Energy;
    const double fm = fhll.Momentum;
    const double fe = fhll.Energy;
    */
    const double um = 0.5*(uhll.positive-uhll.negative);
    const double ue = 0.5*(uhll.positive+uhll.negative);
    const double fm = 0.5*(fhll.positive-fhll.negative);
    const double fe = 0.5*(fhll.positive+fhll.negative);
    const double temp = abs((fe+um+ue+fm)/(fe+um-ue-fm));
    if(temp<thres)
      return -sqrt((fe+0.5*(fm+ue))/(fe+fm+ue+um));
    else if(temp>1/thres)
      return sqrt((-fe+0.5*(fm+ue))/(-fe+fm+ue-um));
    const double aux1 = ue*sqrt(pow(1+fm/ue,2)-4*(fe/ue)*(um/ue));
    const double aux2 = fm+ue-2*um+aux1;
    const double aux3 = fm+ue+2*um+aux1;
    return 2*um/sqrt(aux2*aux3);
  }

  /*! \brief Criterion for when to use the linear approximation
    \param left Primitive variables on the left side
    \param right Primitive variables on the right side
    \param eos Equation of state
    \param thres Threshold
    \return True if left and right are close enough, false otherwise
   */
  bool linear_valid(const Primitive& left,
		    const Primitive& right,
		    const EquationOfState& eos,
		    double thres = 1e-6)
  {
    return (fabs(left.Pressure-right.Pressure)<thres*(left.Pressure+right.Pressure)) &&
      ((eos.dp2ba(left.Density, left.Pressure)*thres>abs(left.Celerity) &&
	eos.dp2ba(right.Density, right.Pressure)*thres>abs(right.Celerity)) ||
       (left.Celerity-right.Celerity)/(abs(left.Celerity)+abs(right.Celerity))<thres);
  }

  /*! \brief Calculates the pressure at the interface
    \param uhll Hll conserved variables
    \param fhll Hll fluxes
    \return Pressure
   */
  double calc_hll_ps(const NewConserved& uhll,
		     const NewConserved& fhll)
  {
    const double um = 0.5*(uhll.positive - uhll.negative);
    const double ue = 0.5*(uhll.positive + uhll.negative);
    const double fm = 0.5*(fhll.positive - fhll.negative);
    const double fe = 0.5*(fhll.positive + fhll.negative);
    return -0.5*(ue-fm-sqrt(pow(ue,2)+pow(fm,2)+2*ue*fm-4*um*fe));
  }
}

RiemannSolution HLL::operator()(Primitive const& left_i, 
				Primitive const& right_i) const
{
  if(fabs(left_i.Pressure)<=0 &&
     fabs(right_i.Pressure)<=0 &&
     right_i.Celerity>left_i.Celerity)
    return RiemannSolution(0,0);
  if(fabs(left_i.Pressure)<=0 &&
     fabs(right_i.Pressure)<=0 &&
     fabs(left_i.Celerity-right_i.Celerity)<=0)
    return RiemannSolution(0,0);
  if(linear_valid(left_i, right_i, eos_))
    return (LinearRS(eos_))(left_i,right_i);
  /*
  const NewConserved new_uhll = CalcUhll(left_i, right_i, eos_);
  const NewConserved new_fhll = CalcFhll(left_i, right_i, eos_);
  */
  const double offset = 0.5*(left_i.Celerity + right_i.Celerity);
  const Primitive left(left_i.Density,
		       left_i.Pressure,
		       celerity_addition(left_i.Celerity,
					 -offset));
  const Primitive right(right_i.Density,
			right_i.Pressure,
			celerity_addition(right_i.Celerity,
					  -offset));
  const std::pair<NewConserved,NewConserved> uhll_fhll = 
    calc_uhll_fhll(left, right, eos_);
  const NewConserved& uhll = uhll_fhll.first;
  const NewConserved& fhll = uhll_fhll.second;
  const double ws = calc_hll_ws(uhll,fhll);
  const double ps = calc_hll_ps(uhll, fhll);

  return RiemannSolution(ps,celerity_addition(ws,offset));
}
