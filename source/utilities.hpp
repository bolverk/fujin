/*! \file utilities.hpp
  \brief Various usefull utility functions
  \author Almog Yalinewich
 */

#ifndef UTILITIES_HPP
#define UTILITIES_HPP 1

#include <complex>
#include <vector>
#include "trans_eqn_solver.hpp"

using std::complex;
using std::string;
using std::vector;

/*! \brief Relativistic velocity addition
  \param v1 First velocity
  \param v2 Second velocity
  \return Velocity
 */
double RelVelAdd(double v1, double v2);

/*! \brief Differential of the velocity addition
  \param b1 First velocity
  \param b2 Second velocity
  \param db1 Differential of the first velocity
  \param db2 Differential of the second velocity
 */
double dRelVelAdd(double b1, double b2, double db1, double db2);

/*! \brief Celerity addition
  \param w1 First celerity
  \param w2 Second celerity
  \return Sum of celerities
 */
double celerity_addition(double w1, double w2);

/*! \brief Differential of the celerity addition
  \param w1 First celerity
  \param w2 Second celerity
  \param dw1 Differential of the first celerity
  \param dw2 Differential of the second celerity
  \return Differential of sum of celerities
 */
double celerity_addition_diff(double w1, double w2,
			      double dw1, double dw2);

/*! \brief Converts velocity to celerity
  \param v Velocity
  \return Celerity
 */
double velocity2celerity(double v);

/*! \brief Converts celerity to velocity
  \param w Celerity
  \return Velocity
 */
double celerity2velocity(double w);

/*! \brief Calculates the Lorentz factor
  \param w Celerity
  \return Loretnz factor
 */
double celerity2lorentz_factor(double w);

/*! \brief Calculates the Loretnz factor
  \param Velocity Velocity
  \result Lorentz factor
 */
double Velocity2LorentzFactor(double Velocity);

/*! \brief Multiplication of a complex number with an integer
  \param a Integer
  \param c Complex number
  \return Product
 */
complex<double> operator*(int a, complex<double> c);

/*! \brief Subtraction of an integer and a complex number
  \param a Integer
  \param c Complex number
  \return Difference
 */
complex<double> operator-(int a, complex<double> c);

/*! \brief Sum of an integer and a complex number
  \param a Integer
  \param c Complex number
  \return Sum
 */
complex<double> operator+(int a, complex<double> c);

/*! \brief Converts an integer into a string
  \param n Number
  \return Number as a string
 */
string int2str(int n);

/*! \brief Check whether a variable contains a nan
  \param x A number
  \return True if x is a NaN, 0 otherwise
 */
bool is_nan(double x);

/*! \brief Calculates the smallest term in an array
  \param v Array
  \return Smallest term
 */
double min(vector<double> const& v);

/*! \brief Calculates the maximum term
  \param v List of numbers
  \return Maximum term
 */
template<class MT> MT max_term(vector<MT> const& v)
{
  MT res = v[0];
  for(size_t i=1, endp=v.size();i<endp;++i)
    res = std::max(res,v[i]);
  return res;
}

/*! \brief Calculates the minium term
  \param v List of numbers
  \return Minimum term
 */
template<class MT> MT min_term(vector<MT> const& v)
{
  MT res = v[0];
  for(size_t i=1, endp=v.size();i<endp;++i)
    res = std::min(res,v[i]);
  return res;
}

/*! \brief Array of uniformly spaced real values
  \param vmin Minimum value
  \param vmax Maximum value
  \param num Number of terms
  \return Vector of uniformly spaced values
 */
vector<double> linspace(double vmin, double vmax, size_t num);

vector<double> logspace(double vmin, double vmax, size_t num,
			double dv_ratio);

//! \brief A scalar function
class ScalarFunction 
{
public:

  /*! \brief Evaluates the function
    \param x Argument
    \return Result
   */
  virtual double Eval(double x) const = 0;

  virtual ~ScalarFunction(void);
};

/*! \brief Applies a scalar function to all members of a vector
  \param v Vector
  \param sf Scalar function
  \return Modified vector
 */
vector<double> apply_to_all_members(vector<double> const& v,
				    ScalarFunction const& sf);

vector<double> operator*(double d, vector<double> const& v);

vector<double> join(vector<double> const& v1,
		    vector<double> const& v2);

/*! \brief Checks if a number is effectively zero withing working precision
  \param x Number
  \return True if x if very close to zero
 */
bool effectively_zero(double x);

//! \brief Lazy list (evaluated on demand)
template<class T> class Index2Member
{
public:

  /*! \brief Returns the length of the list
    \return Length of the list
   */
  virtual size_t getLength(void) const = 0;

  /*! \brief Evaluates list member
    \param i Index
    \return i'th list member
   */
  virtual T operator()(size_t i) const = 0;

  virtual ~Index2Member(void) {}
};

//! \brief Converts a vector to Index2Member
template<class T> class Echo: public Index2Member<T>
{
public:

  /*! \brief Class constructor
    \param v stl vector
   */
  Echo(const vector<T>& v):
    v_(v) {}

  size_t getLength(void) const
  {
    return v_.size();
  }

  T operator()(size_t i) const
  {
    return v_[i];
  }

private:
  //! \brief Reference to stl vector
  const vector<T>& v_;
};

/*! \brief Generates a vector
  \param i2m Lazy list
  \return stl vector
 */
template<class T> vector<T> serial_generate(const Index2Member<T>& i2m)
{
  vector<T> res(i2m.getLength());
  for(size_t i=0, endp=res.size();i<endp;++i)
    res[i] = i2m(i);
  return res;
}

/*! \brief Sum of all members
  \param i2m Lazy list
  \return Sum of all members
 */
template<class T> T sum_all(const Index2Member<T>& i2m)
{
  T res = i2m(0);
  for(size_t i=1;i<i2m.getLength();++i)
    res += i2m(i);
  return res;
}

/*! \brief Addition operator
  \param v1 Left argument
  \param v2 Right argument
  \return Sum
 */
template<class T> vector<T> operator+(const vector<T>& v1,
				      const vector<T>& v2)
{
  assert(v1.size()==v2.size());
  vector<T> res(v1.size());
  for(size_t i=0;i<v2.size();++i)
    res[i] = v1[i] + v2[i];
  return res;
}

/*! \brief Subtraction operator
  \param v1 Left argument
  \param v2 Right argument
  \return Difference
 */
template<class T> vector<T> operator-(const vector<T>& v1,
				      const vector<T>& v2)
{
  assert(v1.size()==v2.size());
  vector<T> res(v1.size());
  for(size_t i=0;i<v1.size();++i)
    res[i] = v1[i] - v2[i];
  return res;
}

/*! \brief Multiplication by a scalar
  \param s Scalar
  \param v List
  \return Product
 */
template<class S, class T> vector<T> operator*
(const S s, const vector<T>& v)
{
  vector<T> res(v.size());
  for(size_t i=0;i<v.size();++i)
    res[i] = s*v[i];
  return res;
}

/*! \brief min Finds minimum member
  \param i2m Lazy list
  \return Minimum term
 */
template<class T> T min(const Index2Member<T>& i2m)
{
  T res = i2m(0);
  for(size_t i=1;i<i2m.getLength();++i)
    res = min(res,i2m(i));
  return res;
}

#endif // UTILITIES_HPP
