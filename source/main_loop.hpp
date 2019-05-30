#ifndef MAIN_LOOP_HPP
#define MAIN_LOOP_HPP 1

#include <boost/ptr_container/ptr_vector.hpp>
#include "srhd_sim.hpp"

//! \brief Base class for the condition for when to terminate the simulation
class TerminationCondition
{
public:
  
  /*! \brief Returns true if the simulation should stop
    \param sim Reference to the simulation
    \return True if simulation should stop
   */
  virtual bool operator()(SRHDSimulation const& sim) const = 0;

  virtual ~TerminationCondition(void);
};

//! \brief Terminates the simulation after a certain number of iterations
class IterationTermination: public TerminationCondition
{
public:

  /*! \brief Class constructor
    \param max_iter Number of iterations
   */
  explicit IterationTermination(int max_iter);

  bool operator()(SRHDSimulation const& sim) const;

private:

  //! \brief Maximum number of iterations
  const int max_iter_;
};

//! \brief Terminates the simulation after a certain amount of virtual time
class SafeTimeTermination: public TerminationCondition
{
public:

  /*! \brief Class constructor
    \param tf Termination time
    \param max_iter Upper limit on number of iterations
   */
  SafeTimeTermination(double tf,
		      int max_iter);

  bool operator()(SRHDSimulation const& sim) const;

private:

  //! \brief Termination time
  const double tf_;

  //! \brief Upper limit on number of iterations
  const int max_iter_;
};

//! \brief Base class for diagnostic function
class DiagnosticFunction
{
public:

  /*! \brief Perform diagnostic
    \param sim Reference to simulation
   */
  virtual void operator()(SRHDSimulation const& sim) const = 0;

  virtual ~DiagnosticFunction(void);
};

//! \brief Writes the time each step
class WriteTime: public DiagnosticFunction
{
public:

  /*! \brief Class constructor
    \param fname File name
   */
  explicit WriteTime(string const& fname);

  void operator()(SRHDSimulation const& sim) const;

private:

  //! \brief File name
  const string fname_;
};

//! \brief Returns the history of the energy 
class TotalEnergyHistory: public DiagnosticFunction
{
public:

  /*! \brief Class constructor
    \param fname File name
   */
  explicit TotalEnergyHistory(string const& fname);

  void operator()(SRHDSimulation const& sim) const;

  ~TotalEnergyHistory(void);

private:

  //! \brief fname_ File name
  const string fname_;

  //! \brief energy_ History of the energy
  mutable vector<double> energy_;
};

//! \brief Container for multiple diagnostics
class MultipleDiagnostics: public DiagnosticFunction
{
public:

  MultipleDiagnostics(void);

  void operator()(SRHDSimulation const& sim) const;

  //! \brief List of diagnostics
  boost::ptr_vector<DiagnosticFunction> diag_list;
};

void main_loop(SRHDSimulation& sim,
	       TerminationCondition const& term_cond,
	       void (SRHDSimulation::*time_advance_method)(void),
	       DiagnosticFunction const& diag_func);

#endif // MAIN_LOOP_HPP
