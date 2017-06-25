#include <cmath>
#include <assert.h>
#include "ideal_gas_hydrodynamics.hpp"
#include "ideal_gas.hpp"
#include "utilities.hpp"
#include "trans_eqn_solver.hpp"
#include "advanced_hydrodynamic_variables.hpp"

using namespace std;

double calc_hugoniot_density(double p,
			     double d0,
			     double p0,
			     double g)
{
  const double a = p/d0 + ((d0 + d0*g + g*p)*p0)/(pow(d0,2)*(-1 + g)) + \
    (g*pow(p0,2))/(pow(d0,2)*pow(-1 + g,2));
  const double b = -(((1 + g)*p)/(-1 + g)) - p0;
  const double c = -((g*pow(p,2))/pow(-1 + g,2)) - (g*p*p0)/(-1 + g);
  return (-b+sqrt(pow(b,2)-4*a*c))/(2*a);
}

double calc_hugoniot_density_dp(double p,
				double d0,
				double p0,
				double g)
{
  const double d = calc_hugoniot_density(p,d0,p0,g);
  return (d0*(d*(-1 + g)*(d + d0 - d*g + d0*g) + 2*d0*g*p) + 
	  (-pow(d,2) + pow(d0,2))*(-1 + g)*g*p0)/
    (-(pow(d0,2)*(-1 + g)*((1 + g)*p + (-1 + g)*p0)) + 
     2*d*(g*p0*((-1 + g)*p + p0) + 
	  d0*(-1 + g)*((-1 + g)*p + p0 + g*p0)));
}

double calc_hugoniot_celerity(double p,
			      double d0,
			      double p0,
			      double g)
{
  assert(g>1 && d0>0 && p0>=0 && p>p0 && "Invalid value of the pressure");

  if(fabs(p-p0)/(p+p0)<1e-6)
    return (p-p0)*sqrt((1+d0*(g-1)/(p0*g))*(g-1))/(d0*(g-1)+g*p0);

  const double d = calc_hugoniot_density(p,d0,p0,g);
  const double numerator = (g-1)*(p-p0)*((g-1)*(d-d0)+p-p0);
  const double denominator = (g*p+d*(g-1))*(d0*(g-1)+g*p0);
  return sqrt(numerator/denominator);
}

double calc_hugoniot_celerity_dp(double p,
				 double d0,
				 double p0,
				 double g)
{
  const double d = calc_hugoniot_density(p,d0,p0,g);
  const double dddp = calc_hugoniot_density_dp(p,d0,p0,g);
  const double w = calc_hugoniot_celerity(p,d0,p0,g);
  const double dlnwdp = (pow(d,2)*pow(-1 + g,2) - d*(-1 + g)*
			 (d0*(-1 + g) - 2*p + 2*p0 - g*p0) + 
			 g*(pow(p,2) - p0*(d0*(-1 + g) + p0)))/
    (2.*(d*(-1 + g) + g*p)*(p - p0)*(d0 + d*(-1 + g) - d0*g + p - p0));
  const double dlnwdd = ((-1 + g)*(d0*(-1 + g) + (-1 + g)*p + p0))/
    (2.*(d*(-1 + g) + g*p)*(d0 + d*(-1 + g) - d0*g + p - p0));
  return w*(dlnwdp+dlnwdd*dddp);
}

double calc_isentrope_density(double p,
			      double d0,
			      double p0,
			      double g)
{
  return d0*pow(p/p0,1.0/g);
}

double calc_isentrope_celerity(double p,
			       double d0,
			       double p0,
			       double g)
{
  const IdealGas eos(g);
  const double ba0 = eos.dp2ba(d0,p0);
  const double ba = eos.dp2ba
    (calc_isentrope_density(p,d0,p0,g),p);
  return sinh((2.0/sqrt(g-1))*(atanh(ba/sqrt(g-1))-atanh(ba0/sqrt(g-1))));
}

double calc_isentrope_celerity_dp(double p,
				  double d0,
				  double p0,
				  double g)
{
  const IdealGas eos(g);
  const double ba0 = eos.dp2ba(d0,p0);
  const double ba = eos.dp2ba
    (calc_isentrope_density(p,d0,p0,g),p);
  return (ba/(g*p))*
    cosh((2.0/sqrt(g-1))*(atanh(ba/sqrt(g-1))-atanh(ba0/sqrt(g-1))));
}

double calc_hydrodynamic_celerity(double p,
				  double d0,
				  double p0,
				  double g)
{
  return p>p0 ? calc_hugoniot_celerity(p,d0,p0,g)
    : calc_isentrope_celerity(p,d0,p0,g);
}

double calc_hydrodynamic_celerity_dp(double p,
				     double d0,
				     double p0,
				     double g)
{
  return p>p0 ? calc_hugoniot_celerity_dp(p,d0,p0,g)
    : calc_isentrope_celerity_dp(p,d0,p0,g);
}

double right_celerity(double p,
		      const Primitive& hv,
		      double g)
{
  const double celerity = hv.Celerity;
  return celerity_addition
    (calc_hydrodynamic_celerity(p,hv.Density,hv.Pressure,g),
     celerity);
}

double right_celerity_diff(double p,
			   const Primitive& hv,
			   double g)
{
  const double w1 = hv.Celerity;
  const double dw1 = 0;
  const double w2 = calc_hydrodynamic_celerity
    (p,hv.Density,hv.Pressure,g);
  const double dw2 = calc_hydrodynamic_celerity_dp
    (p,hv.Density,hv.Pressure,g);
  return celerity_addition_diff(w1,w2,dw1,dw2);
}

double left_celerity(double p,
		     Primitive const& hv,
		     double g)
{
  const double celerity = hv.Celerity;
  return celerity_addition
    (-calc_hydrodynamic_celerity(p,hv.Density,hv.Pressure,g),
     celerity);
}

double left_celerity_diff(double p,
			  const Primitive& hv,
			  double g)
{
  const double w1 = hv.Celerity;
  const double dw1 = 0;
  const double w2 = -calc_hydrodynamic_celerity
    (p,hv.Density,hv.Pressure,g);
  const double dw2 = -calc_hydrodynamic_celerity_dp
    (p,hv.Density,hv.Pressure,g);
  return celerity_addition_diff(w1,w2,dw1,dw2);
}

double eval_trans_eqn(double p,
		      Primitive const& left,
		      Primitive const& right,
		      double g)
{
  return right_celerity(p,right,g)-
    left_celerity(p,left,g);
}

double eval_trans_eqn_diff(double p,
			   Primitive const& left,
			   Primitive const& right,
			   double g)
{
  return right_celerity_diff(p,right,g)-
    left_celerity_diff(p,left,g);
}

namespace {

  //! \brief Transcendental equation for the celerity
  class CelerityTranscendentalEquation: public SVDifferentiable
  {
  public:

    /*! \brief Class constructor
      \param left Hydrodynamic variables on the left side of the interface
      \param right Hydrodynamic variables on the right side of the interface
      \param g Adiabatic index
     */
    CelerityTranscendentalEquation(const Primitive& left,
				   const Primitive& right,
				   double g):
      left_(left), right_(right), g_(g) {}

    double operator()(double p) const
    {
      return eval_trans_eqn(p,left_,right_,g_);
    }

    double diff(double p) const
    {
      return eval_trans_eqn_diff(p,left_,right_,g_);
    }

  private:

    //! \brief Hydrodynamic variables on the left side of the interface
    const Primitive left_;

    //! \brief Hydrodynamic variables on the right side of the interface
    const Primitive right_;

    //! \brief Adiabatic index 
    const double g_;
  };
}

double calc_isentrope_pressure(double w,
			       double d0,
			       double p0,
			       double g)
{
  const double ba0 = IdealGas(g).dp2ba(d0,p0);
  const double ba = sqrt(g-1)*
    tanh(0.5*sqrt(g-1)*asinh(w)+atanh(ba0/sqrt(g-1)));
  return p0*pow((g*p0/d0)*(pow(ba,-2)-1.0/(g-1)),-g/(g-1));
}

double overestimate_shock_pressure(double w,
				   double d0,
				   double p0,
				   double g)
{
  const double ba0 = IdealGas(g).dp2ba(d0,p0);
  return p0*(2*pow(ba0,2)+w*g*(w*g+sqrt(4*pow(ba0,2)+pow(w*g,2))))/
    (2*pow(ba0,2));
	     
}

double two_shock_case(const Primitive& left,
		      const Primitive& right,
		      double g)
{
  const CelerityTranscendentalEquation cte(left,right,g);
  const double pmin = 0.5*min(left.Pressure,
			      right.Pressure);
  const double left_celerity = left.Celerity;
  const double right_celerity = right.Celerity;
  const double dw = celerity_addition
    (left_celerity,-right_celerity);
  const double pmax = 2*max
    (overestimate_shock_pressure(dw,
				 left.Density,
				 left.Pressure,g),
     overestimate_shock_pressure(dw,
				 right.Density,
				 right.Pressure,g));
  const SVRelStep sc(1e-3);
  return NRSafe(cte,pmin,pmax,sc);
}

double two_rarefaction_case(const Primitive& left,
			    const Primitive& right,
			    double g)
{
  const CelerityTranscendentalEquation cte(left,right,g);
  const double pmax = 2*max(left.Pressure,
			    right.Pressure);
  const double left_celerity = left.Celerity;
  const double right_celerity = right.Celerity;
  const double dw = celerity_addition
    (left_celerity,-right_celerity);
  const double pmin = 0.5*min
    (calc_isentrope_pressure(dw,
			     right.Density,
			     right.Pressure,g),
     calc_isentrope_pressure(dw,
			     left.Density,
			     left.Pressure,g));
  const SVRelStep sc(1e-3);
  return NRSafe(cte,pmin,pmax,sc);
}

double shock_rarefaction_case(const Primitive& left,
			      const Primitive& right,
			      double g)
{
  const CelerityTranscendentalEquation cte(left,right,g);
  const SVRelStep sc(1e-3);
  const double pmin = min(left.Pressure, right.Pressure)/2;
  const double pmax = max(left.Pressure, right.Pressure)*2;
  return NRSafe(cte,pmin,pmax,sc);
}

namespace {

  /*! \brief Checks whether the acoustic approximation is applicable
    \param left Primitive on the left side
    \param right Primitive on the right
    \param g Adiabatic index
    \param thres Threshold
    \return True if valid, false otherwise
   */
  bool acoustic_approximation_is_valid(const Primitive& left,
				       const Primitive& right,
				       double g,
				       double thres = 1e-3)
  {
    const IdealGas eos(g);
    const double left_sound_speed = eos.dp2ba(left.Density,
					      left.Pressure);
    const double right_sound_speed = eos.dp2ba(right.Density,
					       right.Pressure);
    const bool cond_1 = abs(left.Pressure-right.Pressure)/
      (left.Pressure+right.Pressure)<thres;
    const bool cond_2 = abs(left.Celerity - right.Celerity)/
      (left_sound_speed+right_sound_speed)<thres;
    return cond_1&&cond_2;
  }

  /*! \brief Acoustic approximation to the Riemann problem
    \param left Primitive variables on the left side
    \param right Primitive variables on the right side
    \param g Adiabatic index
    \return Pressure
   */
  double acoustic_approximation(const Primitive& left,
				const Primitive& right,
				double g)
  {
    const IdealGas eos(g);
    const double bal = eos.dp2ba(left.Density,
				 left.Pressure);
    const double bar = eos.dp2ba(right.Density,
				 right.Pressure);
    const double bl = celerity2velocity(left.Celerity);
    const double br = celerity2velocity(right.Celerity);
    const double pl = left.Pressure;
    const double pr = right.Pressure;
    return g*(bl-br+(1-pow(bl,2))*bal/g+(1-pow(br,2))*bar/g)/
      ((1-pow(bl,2))*bal/pl+(1-pow(br,2))*bar/pr);
  }
}

double solve_trans_eqn(Primitive const& left,
		       Primitive const& right,
		       double g)
{
  if(acoustic_approximation_is_valid(left, right, g))
    return acoustic_approximation(left, right, g);

  const double wrpl = right_celerity
    (left.Pressure, right, g);
  const double wlpr = left_celerity
    (right.Pressure, left, g);
  const double right_celerity = right.Celerity;
  const double left_celerity = left.Celerity;
  if(max(wrpl,right_celerity)<min(wlpr,left_celerity))
    return two_shock_case(left,right,g);
  else if(min(wrpl,right_celerity)>max(wlpr,left_celerity))
    return two_rarefaction_case(left,right,g);
  else
    return shock_rarefaction_case(left,right,g);
}
