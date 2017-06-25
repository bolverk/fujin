#include <iostream>
#include <cmath>
#include "collela.hpp"

using namespace std;

Collela::Collela(double ref, double alpha, 
		 double l, double x0):
  ref_(ref), a_(alpha), l_(l), x0_(x0) {}

double Collela::operator()(double x) const
{
  const double xc = x - x0_;
  double f = 0;
  if(abs(xc)<l_)
    f = pow(pow(xc/l_,2)-1,4);
  return ref_*(1+a_*f);
}
