/*
  Tests the function 
  HydrodynamicVelocity from igrs
 */

#include <iostream>
#include <fstream>
#include <vector>
#include "igrs.hpp"

using namespace std;

int main()
{
  const double p0 = 1.;
  const double pmin = 0.05;
  const double pmax = 20.;
  const double g = 5./3.;
  const double bmin = HydrodynamicVelocity(pmin,p0,g);
  const double bmax = HydrodynamicVelocity(pmax,p0,g);

  // Write data to file
  ofstream f("res.txt");
  f << bmin << endl;
  f << bmax << endl;
  f.close();

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
