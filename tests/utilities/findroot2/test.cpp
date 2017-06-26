/*
  Tests the function 
  NRSafe
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
  MyFunc func;
  ofstream f;
  f.open("res.txt");
  const SVRelStep sc(1e-9);
  f << NRSafe(func,0.2,  2., sc);
  f.close();

  ofstream("test_terminated_normally.res").close();
 return 0;
}
