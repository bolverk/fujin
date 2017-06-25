/* \brief Implementation of the functions in riemann_solver.hpp
   \file riemann_solver.cpp
   \author Almog Yalinewich
*/

#include "riemann_solver.hpp"
#include "utilities.hpp"

using namespace std;

RiemannSolution::RiemannSolution(void):
  Pressure(0), Celerity(0) {}

RiemannSolution::RiemannSolution(double p, double w):
  Pressure(p), Celerity(w) {}

RiemannSolver::~RiemannSolver(void) {}
