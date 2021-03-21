/*! \file trigger.hpp
\brief Trigger for diagnostics
\author Almog Yalinewich
 */

#ifndef TRIGGER_HPP
#define TRIGGER_HPP 1

#include "srhd_sim.hpp"

//! \brief Trigger for diagnostics
#if SCAFFOLDING != 1
template<class CE, class CP>
#endif // SCAFFOLDING
class Trigger
{
public:

  /*! \brief returns true when some condition is satisfied
    \param sim Simulation
    \return True if condition is met, false otherwise
   */
#if SCAFFOLDING == 1
  virtual bool operator()(const SRHDSimulation& sim) = 0;
#else 
  virtual bool operator()(const SRHDSimulation<CE, CP>& sim) = 0;
#endif // SCAFFOLDING
  
  /*! \brief Returns the number of times trigger returned true
    \return Number of times trigger returned true
   */
  virtual int getCount(void) const = 0;

  //! \brief Class destructor
  virtual ~Trigger(void);
};

#endif // TRIGGER_HPP
