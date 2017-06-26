/*
  Tests the function 
  FindRoot from igrs
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "trans_eqn_solver.hpp"

using namespace std;

namespace {
  class MyFunc: public SVDifferentiable
  {
  public:
    MyFunc(void) {}

    double operator()(double x) const
    {
      return pow(x,3)-1;
    }
    double diff(double x) const
    {
      return 3*pow(x,2);
    }
  };
}

int main()
{
  const MyFunc func;
  const double xg = 2.;
  const SVRelStep sc(1e-12);
  const double res = NewtonRaphson(func, xg, sc);

  ofstream f("res.txt");
  f << res << endl;
  f.close();

  ofstream("test_terminated_normally.res").close();
 return 0;
}
