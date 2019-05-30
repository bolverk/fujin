#include "aimm_recovery.hpp"
#include <cmath>
#include <cassert>
#include <iostream>

namespace {

  std::pair<double,double> eval_func_deriv(const NewConserved& c,
					   double g,
					   double p)
  {
    const double tau2d = 0.5*(c.positive + c.negative) - 1;
    const double p2d = p/c.mass;
    const double Ws = (tau2d+p2d+1)/sqrt((c.positive+p2d)*(c.negative+p2d));
    const double ds = c.mass/Ws;
    const double es = (tau2d+1-Ws+p2d*(1-pow(Ws,2)))/Ws;
    const double vs2 = 1. - 1./pow(Ws,2);
    const double cs2 = (g-1)*g*es/(1+g*es);
    return std::pair<double,double>((g-1)*ds*es-p,
				    vs2*cs2-1);
  }
}

double calc_pressure(const NewConserved& c,
		     double g)
{
  double prev = fmax(0,-c.mass*fmin(c.positive,c.negative));
  std::pair<double, double> f_df = eval_func_deriv(c,g,prev);
  double res = prev - f_df.first / f_df.second;
  int counter = 0;
  while((fabs(f_df.first) > res*1e-9 && fabs(f_df.first)>1e-14) ||
	counter<4){
		if(counter>18){
			std::cout << c.mass << std::endl;
			std::cout << c.positive << std::endl;
			std::cout << c.negative << std::endl;
		}
		assert(++counter<20 && "too many iterations");
//    assert(++counter<20 && "too many iterations");
    //prev = res;
    f_df = eval_func_deriv(c,g,res);
    res -= f_df.first/f_df.second;
  }
  return res;
}
