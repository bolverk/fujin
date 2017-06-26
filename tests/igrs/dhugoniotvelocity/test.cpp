/*
  Tests the function 
  dHugoniotVelocity from igrs
 */

#include <iostream>
#include <fstream>
#include <vector>
#include "igrs.hpp"

using namespace std;

int main()
{
  // Configurations
  const double p0 = 1.;
  const double pm = 2;
  const double pl = 2 - 0.01;
  const double ph = 2 + 0.01;
  const double g = 5./3.;
  
  // Calculate values
  const double bl = HugoniotVelocity(pl,p0,g);
  const double bh = HugoniotVelocity(ph,p0,g);
  const double dba = dHugoniotVelocity(pm,p0,g);
  const double dbn = (bh-bl)/(ph-pl);

  // Write data to file
  ofstream f("res.txt");
  f << dba << endl;
  f << dbn << endl;
  f.close();

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
