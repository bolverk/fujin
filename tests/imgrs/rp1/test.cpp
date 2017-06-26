/*
  Checks how the solver handles a Riemann problem that caused the simulation to crash
 */

#include <iostream>
#include <fstream>
#include "imgrs.hpp"
#include "ideal_gas.hpp"
#include "advanced_hydrodynamic_variables.hpp"
#include "utilities.hpp"

using namespace std;

int main()
{
  const IdealGas eos(4./3.);
  const IdealGasRiemannSolver rs(eos.getAdiabaticIndex());
  const Primitive left(0.0001305, 0.0243256, 0);
  const Primitive right(0.000151612, 0.0158653, velocity2celerity(0.670679));
  const RiemannSolution res = rs(left,right);

  ofstream f("res.txt");
  f << res.Pressure << endl;
  f << celerity2velocity(res.Celerity) << endl;
  f.close();

  ofstream("test_terminated_normally.res").close();
 return 0;
}
