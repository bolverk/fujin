#include <boost/foreach.hpp>
#include "main_loop.hpp"
#include "universal_error.hpp"
#include "diagnostics.hpp"

TerminationCondition::~TerminationCondition(void) {}

IterationTermination::IterationTermination(int max_iter):
  max_iter_(max_iter) {}

bool IterationTermination::operator()(SRHDSimulation const& sim) const
{
  return sim.GetCycle()>max_iter_;
}

SafeTimeTermination::SafeTimeTermination
(double tf, int max_iter):
  tf_(tf), max_iter_(max_iter) {}

bool SafeTimeTermination::operator()(SRHDSimulation const& sim) const
{
  if(max_iter_>0)
    assert(sim.GetCycle()<max_iter_);
  return sim.GetTime()>tf_;
}

DiagnosticFunction::~DiagnosticFunction(void) {}

WriteTime::WriteTime(string const& fname):
  fname_(fname) {}

void WriteTime::operator()(SRHDSimulation const& sim) const
{
  write_number(sim.GetTime(),fname_.c_str());
}

TotalEnergyHistory::TotalEnergyHistory(string const& fname):
  fname_(fname), energy_() {}

void TotalEnergyHistory::operator()(SRHDSimulation const& sim) const
{
  energy_.push_back(TotalEnergy(sim));
}

TotalEnergyHistory::~TotalEnergyHistory(void)
{
  write_to_file(energy_,fname_.c_str());
}

MultipleDiagnostics::MultipleDiagnostics(void):
  diag_list() {}

void MultipleDiagnostics::operator()(SRHDSimulation const& sim) const
{
  for(size_t i=0;i<diag_list.size();++i)
    diag_list[i](sim);
}

void main_loop(SRHDSimulation& sim,
	       TerminationCondition const& term_cond,
	       void (SRHDSimulation::*time_advance_method)(void),
	       DiagnosticFunction const& diag_func)
{
  while(!term_cond(sim)){
    (sim.*time_advance_method)();

    diag_func(sim);
  }
}
