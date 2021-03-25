#ifndef MAIN_LOOP_HPP
#define MAIN_LOOP_HPP 1

#include <boost/ptr_container/ptr_vector.hpp>
#include "srhd_sim.hpp"
#include "diagnostics.hpp"

//! \brief Base class for the condition for when to terminate the simulation
template<class CE, class CP> class TerminationCondition
{
public:
  
  /*! \brief Returns true if the simulation should stop
    \param sim Reference to the simulation
    \return True if simulation should stop
  */
  virtual bool operator()(const SRHDSimulation<CE, CP>& sim) const = 0;

  virtual ~TerminationCondition(void) {}
};

//! \brief Terminates the simulation after a certain number of iterations
template<class CE, class CP> class IterationTermination: public TerminationCondition<CE, CP>
{
public:

  /*! \brief Class constructor
    \param max_iter Number of iterations
  */
  explicit IterationTermination(int max_iter):
    max_iter_(max_iter) {}

  bool operator()(const SRHDSimulation<CE, CP>& sim) const
  {
    return sim.getCycle()>max_iter_;
  }

private:

  //! \brief Maximum number of iterations
  const int max_iter_;
};

//! \brief Terminates the simulation after a certain amount of virtual time
template<class CE, class CP> class SafeTimeTermination: public TerminationCondition<CE, CP>
{
public:

  /*! \brief Class constructor
    \param tf Termination time
    \param max_iter Upper limit on number of iterations
  */
  SafeTimeTermination(double tf,
		      int max_iter):
    tf_(tf), max_iter_(max_iter) {}

  bool operator()(const SRHDSimulation<CE, CP>& sim) const
  {
    if(max_iter_>0)
      assert(sim.getCycle()<max_iter_);
    return sim.getTime()>tf_;
  }

private:

  //! \brief Termination time
  const double tf_;

  //! \brief Upper limit on number of iterations
  const int max_iter_;
};

//! \brief Base class for diagnostic function
template<class CE, class CP> class DiagnosticFunction
{
public:

  /*! \brief Perform diagnostic
    \param sim Reference to simulation
  */
  virtual void operator()(const SRHDSimulation<CE, CP>& sim) const = 0;

  virtual ~DiagnosticFunction(void) {}
};

//! \brief Writes the time each step
template<class CE, class CP> class WriteTime: public DiagnosticFunction<CE, CP>
{
public:

  /*! \brief Class constructor
    \param fname File name
  */
  explicit WriteTime(const string& fname):
    fname_(fname) {}

  void operator()(const SRHDSimulation<CE, CP>& sim) const
  {
    write_number(sim.getTime(),fname_.c_str());
  }

private:

  //! \brief File name
  const string fname_;
};

//! \brief Returns the history of the energy
template<class CE, class CP> class TotalEnergyHistory: public DiagnosticFunction<CE, CP>
{
public:

  /*! \brief Class constructor
    \param fname File name
  */
  explicit TotalEnergyHistory(string const& fname):
    fname_(fname), energy_() {}

  void operator()(const SRHDSimulation<CE, CP>& sim) const
  {
    energy_.push_back(TotalEnergy(sim));
  }

  ~TotalEnergyHistory(void)
  {
    write_to_file(energy_,fname_.c_str());
  }

private:

  //! \brief fname_ File name
  const string fname_;

  //! \brief energy_ History of the energy
  mutable vector<double> energy_;
};

//! \brief Container for multiple diagnostics
template<class CE, class CP> class MultipleDiagnostics: public DiagnosticFunction<CE, CP>
{
public:

  MultipleDiagnostics(void);

  void operator()(SRHDSimulation<CE, CP> const& sim) const;

  //! \brief List of diagnostics
  boost::ptr_vector<DiagnosticFunction<CE, CP> > diag_list;
};

template <class CE, class CP>
void main_loop
(SRHDSimulation<CE, CP>& sim,
 TerminationCondition<CE, CP> const& term_cond,
 void (SRHDSimulation<CE, CP>::*time_advance_method)(void),
 DiagnosticFunction<CE, CP> const& diag_func)
{
  while(!term_cond(sim)){
    (sim.*time_advance_method)();
    diag_func(sim);
  }
}

#endif // MAIN_LOOP_HPP
