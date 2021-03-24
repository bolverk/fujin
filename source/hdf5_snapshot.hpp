/*! \file hdf5_snapshot.hpp
  \brief Writes hydrodynamic snapshot to hdf5 files
  \author Almog Yalinewich
 */

#ifndef HDF5_SNAPSHOT_HPP
#define HDF5_SNAPSHOT_HPP 1

#include "srhd_sim.hpp"


#if SCAFFOLDING != 1
#include <H5Cpp.h>
#include "hdf5_utils.hpp"
#include "diagnostics.hpp"

using namespace H5;

  //! \brief Calculates the average of every consecutive pair of terms
  class MidValues: public Index2Member<double>
  {
  public:

    /*! \brief Class constructor
      \param i2m Source list
     */
    explicit MidValues(const Index2Member<double>& i2m):
      i2m_(i2m) {}

    size_t getLength(void) const
    {
      return i2m_.getLength() - 1;
    }

    double operator()(size_t i) const
    {
      return 0.5*(i2m_(i)+i2m_(i+1));
    }
  private:

    //! \brief Source list
    const Index2Member<double>& i2m_;
  };
#endif // SCAFFOLDING

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
 const string& fname)
 #if SCAFFOLDING == 1
 ;
 #else
 {
   (HDF5Shortcut(fname))
    ("position",
     serial_generate(MidValues(Echo<double>(sim.getHydroSnapshot().edges))))
	("edges",
     serial_generate(Echo<double>(sim.getHydroSnapshot().edges)))
    ("density",
     serial_generate(PrimitivePropertyGetter<CE, CP>(sim,0)))
    ("pressure",
     serial_generate(PrimitivePropertyGetter<CE, CP>(sim,1)))
    ("celerity",
     serial_generate(PrimitivePropertyGetter<CE, CP>(sim,2)))
    ("time",
     vector<double>(1,sim.getTime()))
    ("cycle",
     vector<double>(1,sim.getCycle()));
 }
 #endif 

/*! \brief Read an hdf5 snapshot
  \param fname File name
  \return Hydrodynamic snapshot
 */
HydroSnapshot read_hdf5_snapshot
(const string& fname);

#endif // HDF5_SNAPSHOT_HPP

