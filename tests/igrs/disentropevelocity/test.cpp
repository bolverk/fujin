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
  const double pm = 0.1;
  const double pl = pm*0.99;
  const double ph = pm*1.01;
  const double g = 5./3.;
  
  // Calculate values
  const double bl = IsentropeVelocity(pl,p0,g);
  const double bh = IsentropeVelocity(ph,p0,g);
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
