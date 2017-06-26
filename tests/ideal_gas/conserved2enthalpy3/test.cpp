/*
  Verifies calculation of enthalpy
 */

#include <iostream>
#include <fstream>
#include <vector>
#include "ideal_gas.hpp"
#include "utilities.hpp"
#include "hydrodynamic_variables.hpp"
#include "advanced_hydrodynamic_variables.hpp"

using namespace std;

int main()
{
  const IdealGas eos(4./3.);

  const Primitive hs(1,1e-6,velocity2celerity(1-1e-6));
  const Conserved cv = Primitive2Conserved(hs,eos);

  const double dh = conserved2enthalpy_diff
    (cv.Energy+cv.Momentum,
     cv.Energy-cv.Momentum,
     eos.getAdiabaticIndex());

  // Write data to file
  ofstream f("res.txt");
  f << eos.dp2e(hs.Density,hs.Pressure) + hs.Pressure-1 << endl;
  f << dh << endl;
  f.close();

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
