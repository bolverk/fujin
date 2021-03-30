/*! \file consecutive_snapshots.hpp
  \brief Writes a series of snapshots to the disk
  \author Almog Yalinewich
*/

//! \brief Writes snapshots at fixed intervals
#include "main_loop.hpp"
#include "trigger.hpp"
#include "filename_pattern.hpp"

template<class CE, class CP>
class ConsecutiveSnapshots: public DiagnosticFunction
<CE, CP>
{
public:

  /*! \brief Class constructor
    \param trigger Trigger function
    \param fnp File name generator
  */
  ConsecutiveSnapshots
   (Trigger<CE, CP>& trigger,
    FileNamePattern& fnp):
     trigger_(trigger), fnp_(fnp) {}    

  /*! \brief perform diagnostic
    \param sim Hydrodynamic simulation
   */
   void operator()(const SRHDSimulation<CE, CP>& sim) const override
  {
    if(trigger_(sim))
      write_hdf5_snapshot(sim, fnp_(trigger_.getCount()));
  }

private:

  //! \brief Snapshot trigger
   Trigger<CE, CP>& trigger_;

  //! \brief Generator for snapshot file name
  FileNamePattern& fnp_;
};
