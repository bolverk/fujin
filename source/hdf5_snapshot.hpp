/*! \file hdf5_snapshot.hpp
  \brief Writes hydrodynamic snapshot to hdf5 files
  \author Almog Yalinewich
 */

#ifndef HDF5_SNAPSHOT_HPP
#define HDF5_SNAPSHOT_HPP 1

#include "srhd_sim.hpp"

/*! \brief Writes simulation data to hdf5
  \param sim Simulation
  \param fname File name
 */
#if SCAFFOLDING != 1
template<class CE, class CP>
#endif // SCAFFOLDING
void write_hdf5_snapshot
#if SCAFFOLDING == 1
(const SRHDSimulation& sim,
#else
 (const SRHDSimulation<CE, CP>& sim,
#endif // SCAFFOLDING
 const string& fname);

/*! \brief Read an hdf5 snapshot
  \param fname File name
  \return Hydrodynamic snapshot
 */
HydroSnapshot read_hdf5_snapshot
(const string& fname);

#endif // HDF5_SNAPSHOT_HPP

