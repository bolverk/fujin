#include "trans_eqn_solver.hpp"
#include <vector>

using std::vector;
using std::size_t;

//! \brief A polynomial
class Polynomial: public SVDifferentiable
{
public:

  /*! \brief Class constructor
    \param coefs List of coefficients
   */
  Polynomial(const vector<double>& coefs);

  double operator()(double x) const;

  double diff(double x) const;

private:

  //! \brief Polynomial coefficients
  const vector<double> coefs_;

  //! \brief Polynomial coefficients of first derivative
  const vector<double> deriv_coefs_;
};
