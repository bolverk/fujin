#ifndef MAIN_LOOP_HPP
#define MAIN_LOOP_HPP 1

#include <boost/ptr_container/ptr_vector.hpp>
#include "srhd_sim.hpp"
#include "diagnostics.hpp"

//! \brief Base class for the condition for when to terminate the simulation
#if SCAFFOLDING == 1

class TerminationCondition
#else
template<class CE, class CP> class TerminationCondition
#endif // SCAFFOLDING
{
public:
  
  /*! \brief Returns true if the simulation should stop
    \param sim Reference to the simulation
    \return True if simulation should stop
   */
#if SCAFFOLDING == 1
  virtual bool operator()(const SRHDSimulation& sim) const = 0;
#else
virtual bool operator()(const SRHDSimulation<CE, CP>& sim) const = 0;
#endif // SCAFFOLDING

#if SCAFFOLDING == 1
  virtual ~TerminationCondition(void);
#else
virtual ~TerminationCondition(void) {}
#endif // SCAFFOLDING
};

//! \brief Terminates the simulation after a certain number of iterations
#if SCAFFOLDING == 1
class IterationTermination: public TerminationCondition
#else
template<class CE, class CP> class IterationTermination: public TerminationCondition<CE, CP>
#endif
{
public:

  /*! \brief Class constructor
    \param max_iter Number of iterations
   */
#if SCAFFOLDING == 1
  explicit IterationTermination(int max_iter);
#else
  explicit IterationTermination(int max_iter):
    max_iter_(max_iter) {}
#endif // SCAFFOLDING

#if SCAFFOLDING == 1
  bool operator()(const SRHDSimulation& sim) const;
#else
  bool operator()(const SRHDSimulation<CE, CP>& sim) const
  {
    return sim.getCycle()>max_iter_;
  }
#endif // SCAFFOLDING

private:

  //! \brief Maximum number of iterations
  const int max_iter_;
};

//! \brief Terminates the simulation after a certain amount of virtual time
#if SCAFFOLDING == 1
class SafeTimeTermination: public TerminationCondition
#else
template<class CE, class CP> class SafeTimeTermination: public TerminationCondition<CE, CP>
#endif // SCAFFOLDING
{
public:

  /*! \brief Class constructor
    \param tf Termination time
    \param max_iter Upper limit on number of iterations
   */
#if SCAFFOLDING == 1
  SafeTimeTermination(double tf,
		      int max_iter);
#else
  SafeTimeTermination(double tf,
		      int max_iter):
    tf_(tf), max_iter_(max_iter) {}
#endif // SCAFFOLDING

#if SCAFFOLDING == 1
  bool operator()(const SRHDSimulation& sim) const;
#else
  bool operator()(const SRHDSimulation<CE, CP>& sim) const
  {
    if(max_iter_>0)
      assert(sim.getCycle()<max_iter_);
    return sim.getTime()>tf_;
  }
#endif // SCAFFOLDING
  

private:

  //! \brief Termination time
  const double tf_;

  //! \brief Upper limit on number of iterations
  const int max_iter_;
};

//! \brief Base class for diagnostic function
#if SCAFFOLDING == 1
class DiagnosticFunction
#else
template<class CE, class CP> class DiagnosticFunction
#endif // SCAFFOLDING
{
public:

  /*! \brief Perform diagnostic
    \param sim Reference to simulation
   */
#if SCAFFOLDING == 1
  virtual void operator()(const SRHDSimulation& sim) const = 0;
#else
virtual void operator()(const SRHDSimulation<CE, CP>& sim) const = 0;
#endif // SCAFFOLDING 

#if SCAFFOLDING == 1
  virtual ~DiagnosticFunction(void);
#else
virtual ~DiagnosticFunction(void) {};
#endif // SCAFFOLDING
};

//! \brief Writes the time each step
#if SCAFFOLDING == 1
class WriteTime: public DiagnosticFunction
#else
template<class CE, class CP> class WriteTime: public DiagnosticFunction<CE, CP>
#endif // SCAFFOLDING
{
public:

  /*! \brief Class constructor
    \param fname File name
   */
#if SCAFFOLDING == 1
  explicit WriteTime(const string& fname);
#else
  explicit WriteTime(const string& fname):
    fname_(fname) {}
#endif // SCAFFOLDING

#if SCAFFOLDING == 1
  void operator()(const SRHDSimulation& sim) const;
#else
  void operator()(const SRHDSimulation<CE, CP>& sim) const
  {
    write_number(sim.getTime(),fname_.c_str());
  }
#endif // SCAFFOLDING

private:

  //! \brief File name
  const string fname_;
};

//! \brief Returns the history of the energy
#if SCAFFOLDING == 1
class TotalEnergyHistory: public DiagnosticFunction
#else
template<class CE, class CP> class TotalEnergyHistory: public DiagnosticFunction<CE, CP>
#endif // SCAFFOLDING
{
public:

  /*! \brief Class constructor
    \param fname File name
   */
#if SCAFFOLDING == 1
    explicit TotalEnergyHistory(string const& fname);
#else
  explicit TotalEnergyHistory(string const& fname):
  fname_(fname), energy_() {}
#endif // SCAFFOLDING

#if SCAFFOLDING == 1
  void operator()(const SRHDSimulation& sim) const;
#else
  void operator()(const SRHDSimulation<CE, CP>& sim) const
  {
    energy_.push_back(TotalEnergy(sim));
  }
#endif // SCAFFOLDING

#if SCAFFOLDING == 1
  ~TotalEnergyHistory(void);
#else
  ~TotalEnergyHistory(void)
  {
    write_to_file(energy_,fname_.c_str());
  }
#endif // SCAFFOLDING

private:

  //! \brief fname_ File name
  const string fname_;

  //! \brief energy_ History of the energy
  mutable vector<double> energy_;
};

//! \brief Container for multiple diagnostics
#if SCAFFOLDING == 1
class MultipleDiagnostics: public DiagnosticFunction
#else
template<class CE, class CP> class MultipleDiagnostics: public DiagnosticFunction<CE, CP>
#endif // SCAFFOLDING
{
public:

  MultipleDiagnostics(void);

#if SCAFFOLDING == 1
  void operator()(SRHDSimulation const& sim) const;
#else
  void operator()(SRHDSimulation<CE, CP> const& sim) const;
#endif // SCAFFOLDING

  //! \brief List of diagnostics
#if SCAFFOLDING == 1
  boost::ptr_vector<DiagnosticFunction> diag_list;
#else
  boost::ptr_vector<DiagnosticFunction<CE, CP> > diag_list;
#endif // SCAFFOLDING
};

#if SCAFFOLDING != 1
template <class CE, class CP>
#endif // SCAFFOLDING
void main_loop
#if SCAFFOLDING == 1
(SRHDSimulation& sim,
 TerminationCondition const& term_cond,
 void (SRHDSimulation::*time_advance_method)(void),
 DiagnosticFunction const& diag_func);
#else
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
#endif // SCAFFOLDING

#endif // MAIN_LOOP_HPP
