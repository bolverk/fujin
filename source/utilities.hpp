/*! \file utilities.hpp
  \brief Various usefull utility functions
  \author Almog Yalinewich
 */

#ifndef UTILITIES_HPP
#define UTILITIES_HPP 1

#include <complex>
#include <vector>
#include <array>
#include "trans_eqn_solver.hpp"
#include <algorithm>

using std::complex;
using std::string;
using std::vector;
using std::transform;
using std::array;

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
//double min(vector<double> const& v);

/*! \brief Array of uniformly spaced real values
  \param vmin Minimum value
  \param vmax Maximum value
  \param num Number of terms
  \return Vector of uniformly spaced values
 */
vector<double> linspace
(const double vmin,
 const double vmax,
 const size_t num);

//! \brief Creats array with uniformly spaced entries
template<size_t N> array<double, N> linspace
(const double vmin, const double vmax)
{
  const double dx = (vmax-vmin)/static_cast<double>(N-1);
  array<double, N> res;
  std::generate(res.begin(),
	   res.end(),
	   [n = 0, &dx, &vmin]() mutable
	   { return vmin+(n++)*dx; });
  return res;
}

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

vector<double> operator*(double d, vector<double> const& v);

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
template<class T, template<class> class C> class Echo: public Index2Member<T>
{
public:

  /*! \brief Class constructor
    \param v stl vector
   */
  explicit Echo(const C<T>& v):
    v_(v) {}

  size_t getLength(void) const override
  {
    return v_.size();
  }

  T operator()(size_t i) const override
  {
    return v_[i];
  }

private:
  //! \brief Reference to stl vector
  const C<T>& v_;
};

template<class T> using simple_vector=vector<T>;

template<class T> void resize_if_necessary(vector<T>& arr, size_t n)
{
  arr.resize(n);
}

template<class T, size_t M> void resize_if_necessary(array<T, M>, size_t /*n*/)
{
  return;
}

/*! \brief Generates a vector
  \param i2m Lazy list
  \return stl vector
 */
template<class T, template<class> class C=simple_vector> C<T>
  serial_generate(const Index2Member<T>& i2m)
{
  //  C<T> res(i2m.getLength());
  C<T> res;
  resize_if_necessary(res, i2m.getLength());
  size_t n = 0;
  std::generate(res.begin(),
	   res.end(),
	   [&n,&i2m]()mutable
	   {return i2m(n++);});
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
  transform(v1.begin(),
	    v1.end(),
	    v2.begin(),
	    res.begin(),
	    [](const T& t1, const T& t2)
	    {return t1+t2;});
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
  transform(v1.begin(),
	    v1.end(),
	    v2.begin(),
	    res.begin(),
	    [](const T& t1, const T& t2)
	    {return t1-t2;});
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
  transform(v.begin(),
	    v.end(),
	    res.begin(),
	    [&s](const T& t)
	    {return s*t;});
  return res;
}

template<class T> T min(const Index2Member<T>& i2m)
{
  T res = i2m(0);
  for(size_t i=1;i<i2m.getLength();++i)
    res = std::min(res, i2m(i));
  return res;
}

#endif // UTILITIES_HPP
