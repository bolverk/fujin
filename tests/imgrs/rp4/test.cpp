/*
  Checks how the solver handles a strong shock and a strong rarefaction
*/

#include <iostream>
#include <fstream>
#include "imgrs.hpp"
#include "ideal_gas.hpp"
#include "advanced_hydrodynamic_variables.hpp"

using namespace std;

int main()
{
  const IdealGas eos(5./3.);
  const IdealGasRiemannSolver rs(eos.getAdiabaticIndex());
  const Primitive left(1.0,1e12,0.0);
  const Primitive right(1.0,1e-2,0.0);

  const RiemannSolution res = rs(left,right);

  ofstream f("res.txt");
  f << res.Pressure << endl;
  f << celerity2velocity(res.Celerity) << endl;
  f.close();

  ofstream("test_terminated_normally.res").close();
 return 0;
}
