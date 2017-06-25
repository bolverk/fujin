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
void write_hdf5_snapshot(const SRHDSimulation& sim,
			 const string& fname);

#endif // HDF5_SNAPSHOT_HPP

