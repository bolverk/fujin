/*! \file utilities.cpp
  \brief Various usefull utility functions
  \author Almog Yalinewich
 */

#include <cmath>
#include <iostream>
#include <sstream>
#include "utilities.hpp"

using namespace std;

double velocity2celerity(double v)
{
  return v/sqrt(1-pow(v,2));
}

double celerity2velocity(double w)
{
  return w/sqrt(1+pow(w,2));
}

double celerity2lorentz_factor(double w)
{
  return sqrt(pow(w,2)+1);
}

ScalarFunction::~ScalarFunction(void) {}

/*! \brief Relativistic velocity addition
  \param v1 First velocity
  \param v2 Second velcoity
  \return Velocity
 */
double RelVelAdd(double v1, double v2)
{
	return (v1+v2)/(1.0+v1*v2);
}

/*! \brief Full differential of the relativistic velocity addition
  \param b1 First velocity
  \param b2 Second velocity
  \param db1 Differential of the first velocity
  \param db2 Differential of the second velocity
  \return Velocity differential
 */
double dRelVelAdd(double b1, double b2, double db1, double db2)
{
  return (db1+db2)/(1+b1*b2)-
    (b1+b2)*(b2*db1+b1*db2)/pow(1+b1*b2,2);
}

double celerity_addition(double w1, double w2)
{
  return w1*sqrt(pow(w2,2)+1)+w2*sqrt(pow(w1,2)+1);
}

double celerity_addition_diff
(double w1, double w2, double dw1, double dw2)
{
  return (sqrt(pow(w1,2)+1)+w1*w2/sqrt(pow(w2,2)+1))*dw2+
    (sqrt(pow(w2,2)+1)+w1*w2/sqrt(pow(w1,2)+1))*dw1;
}

double Velocity2LorentzFactor(double Velocity)
{
  //  return 1.0/sqrt(1.0-pow(Velocity,2));
  return 1./sqrt((1.0-Velocity)*(1.0+Velocity));
}

complex<double> operator*(int a, complex<double> c)
{
  return (a*1.0)*c;
}

complex<double> operator-(int a, complex<double> c)
{
  return (a*1.0)-c;
}

complex<double> operator+(int a, complex<double> c)
{
  return (a*1.0)+c;
}

bool is_nan(double x)
{
	return x!=x;
}

/*! \brief Uniformly spaced grid
  \param vmin Lower bound
  \param vmax Upper bound
  \param num Number of terms
  \return Equally spaced grid
 */
vector<double> linspace(const double vmin,
			const double vmax,
			const size_t num)
{
  const double dx = (vmax-vmin)/static_cast<double>(num-1);
  vector<double> res(num);
  generate(res.begin(),
	   res.end(),
	   [n = 0, &dx, &vmin]() mutable
	   { return vmin+(n++)*dx; });
  return res;
}

/*! \brief Logarithmically spaced grid
  \param vmin Minimum value
  \param vmax Maximum value
  \param num Number of grid points
  \param q ration between consecutive steps
  \return Logarithmically spaced grid
 */
/*
vector<double> logspace(double vmin, double vmax, 
			size_t num, double q)
{
  const double dv0 = (vmax-vmin)/((pow(q,num)-1)/(q-1));
  vector<double> res(num,0);
  res[0] = vmin;
  for(size_t i=1;i<num;++i){
    res[i] = res[i-1] + pow(dv0,static_cast<double>(i));
  }
  return res;
}
*/

string int2str(int n)
{
  stringstream ss;
  ss << n;
  return ss.str();
}

vector<double> apply_to_all_members(vector<double> const& v,
				    ScalarFunction const& sf)
{
  vector<double> res(v.size(),0);
  for(size_t i=0;i<res.size();++i)
    res[i] = sf.Eval(v[i]);
  return res;
}

/*! \brief Product by a scalar
  \param d Scalar
  \param v List
  \return List where each term is multiplied by d
 */
vector<double> operator*(double d, vector<double> const& v)
{
  vector<double> res(v.size(),0);
  for(size_t i=0;i<v.size();++i)
    res[i] = d*v[i];
  return res;
}

/*! \brief Concatenate two vectors
  \param v1 first vector
  \param v2 second vector
  \return Joint vectors
 */
vector<double> join(vector<double> const& v1,
		    vector<double> const& v2)
{
  vector<double> res;
  for(size_t i=0;i<v1.size();++i)
    res.push_back(v1[i]);
  for(size_t i=0;i<v2.size();++i)
    res.push_back(v2[i]);
  return res;
}

bool effectively_zero(double x)
{
  return abs(x)<1e-14;
}
