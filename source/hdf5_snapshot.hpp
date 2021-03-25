/*! \file hdf5_snapshot.hpp
  \brief Writes hydrodynamic snapshot to hdf5 files
  \author Almog Yalinewich
*/

#ifndef HDF5_SNAPSHOT_HPP
#define HDF5_SNAPSHOT_HPP 1

#include "srhd_sim.hpp"
#include <H5Cpp.h>
#include "hdf5_utils.hpp"
#include "diagnostics.hpp"

//using namespace H5;

//! \brief Calculates the average of every consecutive pair of terms
template<class T> class MidValues: public Index2Member<T>
{
public:

  /*! \brief Class constructor
    \param i2m Source list
  */
  explicit MidValues(const Index2Member<T>& i2m):
    i2m_(i2m) {}

  size_t getLength(void) const
  {
    return i2m_.getLength() - 1;
  }

  double operator()(size_t i) const
  {
    return 0.5*(i2m_(i)+i2m_(i+1));
  }

  ~MidValues(void) {}
private:

  //! \brief Source list
  const Index2Member<T>& i2m_;
};

/*! \brief Writes simulation data to hdf5
  \param sim Simulation
  \param fname File name
*/
template<class CE, class CP>
void write_hdf5_snapshot
(const SRHDSimulation<CE, CP>& sim,
 const string& fname)
{
  (HDF5Shortcut(fname))
    ("position",
     serial_generate(MidValues<double>(Echo<double>(sim.getHydroSnapshot().edges))))
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

/*! \brief Read an hdf5 snapshot
  \param fname File name
  \return Hydrodynamic snapshot
*/
HydroSnapshot read_hdf5_snapshot
(const string& fname);

#endif // HDF5_SNAPSHOT_HPP

