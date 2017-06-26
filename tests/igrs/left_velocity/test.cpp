/*
  Tests the functions 
  Right- and LeftVelocity from igrs
 */

#include <iostream>
#include <fstream>
#include <vector>
#include "igrs.hpp"

using namespace std;

int main()
{

  // Configurations
  const double pl = 2.;
  const double p = 1.;
  const double bl = 0.5;
  const double g = 4./3.;

  // Main process
  const double b = LeftVelocity(bl, p, pl ,g);

  // Write data to file
  ofstream f("res.txt");
  f << b << "\n";
  f.close();

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
