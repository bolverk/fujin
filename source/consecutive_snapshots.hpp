/*! \file consecutive_snapshots.hpp
  \brief Writes a series of snapshots to the disk
  \author Almog Yalinewich
*/

//! \brief Writes snapshots at fixed intervals
#include "main_loop.hpp"
#include "trigger.hpp"
#include "filename_pattern.hpp"

#if SCAFFOLDING != 1
template<class CE, class CP>
#endif // SCAFFOLDING
class ConsecutiveSnapshots: public DiagnosticFunction
#if SCAFFOLDING != 1
<CE, CP>
#endif // SCAFFOLDING
{
public:

  /*! \brief Class constructor
    \param trigger Trigger function
    \param fnp File name generator
  */
  ConsecutiveSnapshots
#if SCAFFOLDING == 1
  (Trigger& trigger,
#else
   (Trigger<CE, CP>& trigger,
#endif // SCAFFOLDING		       
    FileNamePattern& fnp);

#if SCAFFOLDING == 1   
  void operator()(const SRHDSimulation& sim) const;
#else
   void operator()(const SRHDSimulation<CE, CP>& sim) const;
#endif // SCAFFOLDING

private:

  //! \brief Snapshot trigger
#if SCAFFOLDING == 1
  Trigger& trigger_;
#else
   Trigger<CE, CP>& trigger_;
#endif // SCAFFOLDING

  //! \brief Generator for snapshot file name
  FileNamePattern& fnp_;
};
