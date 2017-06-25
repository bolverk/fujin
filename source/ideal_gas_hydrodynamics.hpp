#ifndef IDEAL_GAS_HYDRODYNAMICS_HPP
#define IDEAL_GAS_HYDRODYNAMICS_HPP 1

#include "hydrodynamic_variables.hpp"

double calc_hugoniot_density(double p,
			     double d0,
			     double p0,
			     double g);

double calc_hugoniot_density_dp(double p,
				double d0,
				double p0,
				double g);

double calc_hugoniot_celerity(double p,
			      double d0,
			      double p0,
			      double g);

double calc_hugoniot_celerity_dp(double p,
				 double d0,
				 double p0,
				 double g);

double calc_isentrope_density(double p,
			      double d0,
			      double p0,
			      double g);

double calc_isentrope_celerity(double p,
			       double d0,
			       double p0,
			       double g);

double calc_isentrope_celerity_dp(double p,
				  double d0,
				  double p0,
				  double g);

double calc_hydrodynamic_celerity(double p,
				  double d0,
				  double p0,
				  double g);

double calc_hydrodynamic_celerity_dp(double p,
				     double d0,
				     double p0,
				     double g);

double right_celerity(double p,
		      const Primitive& hv,
		      double g);

double right_celerity_diff(double p,
			   const Primitive& hv,
			   double g);

double left_celerity(double p,
		     const Primitive& hv,
		     double g);

double left_celerity_diff(double p,
			  const Primitive& hv,
			  double g);

double eval_trans_eqn(double p,
		      const Primitive& left,
		      const Primitive& right,
		      double g);

double eval_trans_eqn_diff(double p,
			   const Primitive& left,
			   const Primitive& right,
			   double g);

double calc_isentrope_pressure(double w,
			       double d0,
			       double p,
			       double g);

double overestimate_shock_pressure(double w,
				   double d0,
				   double p0,
				   double g);

double two_shock_case(const Primitive& left,
		      const Primitive& right,
		      double g);

double two_rarefaction_case(const Primitive& left,
			    const Primitive& right,
			    double g);

double shock_rarefaction_case(const Primitive& left,
			      const Primitive& right,
			      double g);

double solve_trans_eqn(Primitive const& left,
		       Primitive const& right,
		       double g);

#endif // IDEAL_GAS_HYDRODYNAMICS_HPP
