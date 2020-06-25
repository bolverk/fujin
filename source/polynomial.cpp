#include "polynomial.hpp"

namespace {

  /*! \brief Evaluates a polynomial
    \param x Argumetn
    \param coefs Polynomial coefficients
    \return Value of the polynomial at x
   */
  double poly_eval(double x, const vector<double>& coefs)
  {
    double res = 0;
    for(vector<double>::const_iterator c=coefs.begin(), endp=coefs.end();
	c!=endp;++c)
      res = x*res + (*c);
    return res;	
  }

  /*! \brief Calculates the coefficients of the first derivative
    \param coefs Coefficients of the original polynomial
    \return Coefficients of the first derivative
   */
  vector<double> calc_deriv_coefs(const vector<double>& coefs)
  {
    vector<double> res(coefs.size()-1);
    for(size_t i=0, endp=res.size();i<endp;++i)
      res[i] = coefs[i]*static_cast<double>(endp-i);
    return res;
  }
}

Polynomial::Polynomial(const vector<double>& coefs):
  coefs_(coefs), deriv_coefs_(calc_deriv_coefs(coefs)) {}

double Polynomial::operator()(double x) const
{
  return poly_eval(x, coefs_);
}

double Polynomial::diff(double x) const
{
  return poly_eval(x, deriv_coefs_);
}
