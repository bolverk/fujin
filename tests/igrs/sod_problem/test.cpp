/*
  Tests the function SolveTranscendentalEquation
  The test is based on the Sod problem,
  where matter moves into a rigid wall
  to adapt this problem to the form of a 
  Riemann problem, we use two regions that
  have the same thermodynamic conditions
  and are moving towards each other at
  the same (absolute) velocity
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "igrs.hpp"

using namespace std;

int main()
{

  // Configurations
  const double p0 = 1.;
  const double absb = 0.5;
  const double g = 4./3.;

  // Main process
  const RiemannSolution rs = RiemannSolve(absb, -absb, p0, p0, g);
  const double p = p0*(2*(g-1)+absb*g*sqrt(pow(absb,2)*pow(2-g,2)+
					   4*(g-1))+
		       pow(absb,2)*(pow(g,2)-2*g+2))/
    (2*(1-pow(absb,2))*(g-1));
  const double p2 = SolveTranscendentalEquation(absb, -absb, p0, p0, g);

  // Write data to file
  ofstream f("res.txt");
  f << rs.Velocity << endl;
  f << rs.Pressure << endl;
  f << p << endl;
  f << p2 << endl;
  f.close();

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
