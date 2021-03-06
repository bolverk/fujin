/*
  Verifies that the HLL riemann solver returns the correct answer for the pure advection problem  
 */

#include <iostream>
#include <fstream>
#include <vector>
#include "hll.hpp"
#include "ideal_gas.hpp"
#include "advanced_hydrodynamic_variables.hpp"
#include "utilities.hpp"

using namespace std;

int main()
{
  const IdealGas eos(4./3.);
  const Primitive left(1,0.2,0.1);
  const Primitive right = left;

  // Main process
  const RiemannSolution res = (HLL(eos))(left, right);

  // Write data to file
  ofstream f("res.txt");
  f << res.Pressure << "\n";
  f << res.Celerity << "\n";
  f.close();

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
