/*
  Checks how the solver handles a Riemann problem that caused the simulation to crash
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
  Primitive left(80.6113, 2653.96, 0.452243);
  Primitive right(117.993, 6004.3, -0.0391386);

  const RiemannSolution res1 = rs(left,right);
  left.Celerity = -left.Celerity;
  right.Celerity = -right.Celerity;
  const RiemannSolution res2 = rs(right,left);

  ofstream f("res.txt");
  f << res1.Pressure << endl;
  f << velocity2celerity(res1.Celerity) << endl;
  f << res2.Pressure << endl;
  f << velocity2celerity(res2.Celerity) << endl;
  f.close();

  ofstream("test_terminated_normally.res").close();
 return 0;
}
