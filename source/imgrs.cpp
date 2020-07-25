//! \file imgrs.cpp
//! \author Almog Yalinewich
//! \brief Riemann solver for an ideal, massive, relativistic gas

#include <iostream>
#include <math.h>
#include <limits>
#include "imgrs.hpp"
#include "utilities.hpp"
#include "universal_error.hpp"
#include "ideal_gas_hydrodynamics.hpp"
#include "advanced_hydrodynamic_variables.hpp"

using namespace std;

IdealGasRiemannSolver::IdealGasRiemannSolver(double ig):
  g(ig) {}

RiemannSolution IdealGasRiemannSolver::operator()
  (Primitive const& left_i, Primitive const& right_i) const
{
  if(abs(left_i.Pressure-right_i.Pressure)<std::numeric_limits<double>::epsilon() &&
     abs(left_i.Celerity-right_i.Celerity)<std::numeric_limits<double>::epsilon())
    return RiemannSolution(0.5*(left_i.Pressure+right_i.Pressure),
			   0.5*(left_i.Celerity+right_i.Celerity));
  const double ps = solve_trans_eqn(left_i,right_i,g);
  const double ws = 0.5*(calc_left_celerity(ps,left_i,g)+
			 calc_right_celerity(ps,right_i,g));
  return RiemannSolution(ps,ws);
}
