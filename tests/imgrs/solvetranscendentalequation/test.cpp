/*
  Tests the function
  SolveTranscendentalEquation
*/

#include <iostream>
#include <fstream>
#include <vector>
#include "imgrs.hpp"
#include "ideal_gas.hpp"
#include "ideal_gas_hydrodynamics.hpp"
#include "advanced_hydrodynamic_variables.hpp"

using namespace std;

int main()
{
  const IdealGas eos(4./3.);
  const Primitive left(1.0,0.1,velocity2celerity(0.5));
  const Primitive right(1.0,0.1,velocity2celerity(0.0));

  const double res = solve_trans_eqn(left, right,
				     eos.getAdiabaticIndex());

  ofstream f("res.txt");
  f << res;
  f.close();

  ofstream("test_terminated_normally.res").close();
  return 0;
}
