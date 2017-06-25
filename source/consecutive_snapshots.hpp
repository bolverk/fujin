/*! \file consecutive_snapshots.hpp
  \brief Writes a series of snapshots to the disk
  \author Almog Yalinewich
*/

//! \brief Writes snapshots at fixed intervals
#include "main_loop.hpp"
#include "trigger.hpp"
#include "filename_pattern.hpp"

class ConsecutiveSnapshots: public DiagnosticFunction
{
public:

  /*! \brief Class constructor
    \param trigger Trigger function
    \param fnp File name generator
  */
  ConsecutiveSnapshots(Trigger& trigger,
		       FileNamePattern& fnp);

  void operator()(SRHDSimulation const& sim) const;

private:

  //! \brief Snapshot trigger
  Trigger& trigger_;

  //! \brief Generator for snapshot file name
  FileNamePattern& fnp_;
};
