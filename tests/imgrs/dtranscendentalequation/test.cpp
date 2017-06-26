/*
  Tests the function 
  dTranscendentalEquation
 */

#include <iostream>
#include <fstream>
#include "imgrs.hpp"
#include "ideal_gas.hpp"
#include "ideal_gas_hydrodynamics.hpp"
#include "advanced_hydrodynamic_variables.hpp"

using namespace std;

int main()
{
  // Configurations
  const IdealGas eos(4./3.);
  const Primitive left(1.0, 0.1, velocity2celerity(0.5));
  const Primitive right(1.0, 0.1, velocity2celerity(0.5));
  const double p = 2;
  const double dp = 0.01;
  
  // Calculate values
  const double pl = p - dp;
  const double ph = p + dp;
  const double rl = eval_trans_eqn(p-dp, left, right, 
				   eos.getAdiabaticIndex());
  const double rh = eval_trans_eqn(p+dp, left, right,
				   eos.getAdiabaticIndex());
  const double dra = eval_trans_eqn_diff(p, left, right,
					 eos.getAdiabaticIndex());
  const double drn = (rh-rl)/(ph-pl);

  // Write data to file
  ofstream f("res.txt");
  f << dra << endl;
  f << drn << endl;
  f.close();

  // Finalise
  ofstream("test_terminated_normally.res").close();
 return 0;
}
