#include "consecutive_snapshots.hpp"
#include "hdf5_snapshot.hpp"

ConsecutiveSnapshots::ConsecutiveSnapshots
(Trigger& trigger, FileNamePattern& fnp):
  trigger_(trigger), fnp_(fnp) {}

void ConsecutiveSnapshots::operator()(const SRHDSimulation& sim) const
{
  if(trigger_(sim))
    write_hdf5_snapshot(sim, fnp_(trigger_.getCount()));
}
