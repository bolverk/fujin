/*
  Tests the function
  HugoniotVelocity from igrs
 */

#include <iostream>
#include <fstream>
#include <vector>
#include "igrs.hpp"

using namespace std;

int main()
{
  const double p0 = 1.;
  const double p = 10;
  const double g = 4./3.;
  const double b = HugoniotVelocity(p,p0,g);

  // Write data to file
  ofstream f("res.txt");
  f << b;
  f.close();

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
