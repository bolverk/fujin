/*
  Tests the function 
  dTranscendentalEquation from igrs
 */

#include <iostream>
#include <fstream>
#include <vector>
#include "igrs.hpp"

using namespace std;

int main()
{
  // Configurations
  const double g = 4./3.; 
  const double pl = 2.;
  const double pr = 1.;
  const double bl = 0.5;
  const double br = 0.1;
  
  const double pmid = 0.5*(pl+pr);
  const double pmin = pmid*0.99;
  const double pmax = pmid*1.01;

  // Main process
  const double dea = dTranscendentalEquation(bl, br, pl, pr, pmid, g);
  const double emin = TranscendentalEquation(bl, br, pl, pr, pmin, g);
  const double emax = TranscendentalEquation(bl, br, pl, pr, pmax, g);
  const double den = (emax-emin)/(pmax-pmin);

  // Write data to file
  ofstream f("res.txt");
  f << dea << endl;
  f << den << endl;
  f.close();

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
