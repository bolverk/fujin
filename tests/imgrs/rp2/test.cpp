/*
  Checks how the solver handles a Riemann problem that produced wrong answer
 */

#include <iostream>
#include <fstream>
#include "imgrs.hpp"
#include "ideal_gas.hpp"
#include "advanced_hydrodynamic_variables.hpp"

using namespace std;

int main()
{
  const IdealGas eos(4./3.);
  const IdealGasRiemannSolver rs(eos.getAdiabaticIndex());
  const Primitive left(1.0,1.0,velocity2celerity(0.9));
  const Primitive right(1.0,10.0,0.0);

  const RiemannSolution res = rs(left,right);

  ofstream f("res.txt");
  f << res.Pressure << "\n";
  f << celerity2velocity(res.Celerity) << "\n";
  f.close();

  ofstream("test_terminated_normally.res").close();
 return 0;
}
