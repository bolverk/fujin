/*
  Tests the function
  IsentropeVelocity from igrs
 */

#include <iostream>
#include <fstream>
#include <vector>
#include "igrs.hpp"

using namespace std;

int main()
{
  const double p0 = 1.;
  const double p = 0.1;
  const double g = 4./3.;
  const double b = IsentropeVelocity(p,p0,g);

  // Write data to file
  ofstream f("res.txt");
  f << b << "\n";
  f.close();

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
