/*! \file trigger.hpp
\brief Trigger for diagnostics
\author Almog Yalinewich
 */

#ifndef TRIGGER_HPP
#define TRIGGER_HPP 1

#include "srhd_sim.hpp"

//! \brief Trigger for diagnostics
template<class CE, class CP>
class Trigger
{
public:

  /*! \brief returns true when some condition is satisfied
    \param sim Simulation
    \return True if condition is met, false otherwise
   */
  virtual bool operator()(const SRHDSimulation<CE, CP>& sim) = 0;
  
  /*! \brief Returns the number of times trigger returned true
    \return Number of times trigger returned true
   */
  virtual int getCount(void) const = 0;

  //! \brief Class destructor
  virtual ~Trigger(void);
};

#endif // TRIGGER_HPP
