#include "trans_eqn_solver.hpp"
#include <algorithm>

using std::size_t;
using std::for_each;

//! \brief A polynomial
template<template<class> class C>
class Polynomial: public SVDifferentiable
{
private:

  //! \brief Polynomial coefficients
  const C<double> coefs_;

  //! \brief Polynomial coefficients of first derivative
  const C<double> deriv_coefs_;

  C<double> calc_deriv_coefs_(const C<double>& coefs) const
  {
    C<double> res;
    resize_if_necessary(res, coefs.size());
    for(size_t i=0, endp=res.size()-1;i<endp-1;++i)
      res[i] = coefs[i]*static_cast<double>(endp-i);
    res.front() = 0;
    return res;
  }

  double poly_eval_(double x, const C<double>& coefs) const
  {
    double res = 0;
    for_each(coefs.begin(),
	     coefs.end(),
	     [&res,&x](double c)
	     {res = res*x+c;});
    return res;
  }
public:

  /*! \brief Class constructor
    \param coefs List of coefficients
   */
  explicit Polynomial(const C<double>& coefs):
    coefs_(coefs),
    deriv_coefs_(calc_deriv_coefs_(coefs)) {}

  double operator()(double x) const
  {
    return poly_eval_(x, coefs_);
  }

  double diff(double x) const
  {
    return poly_eval_(x, deriv_coefs_);
  }
};
