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
template<template<class> class CE, template<class> class CP>
void write_hdf5_snapshot
(const SRHDSimulation<CE, CP>& sim,
 const string& fname)
{
  (HDF5Shortcut(fname))
    ("position",
     serial_generate(MidValues<double>(Echo<double, CE>
				       (sim.getHydroSnapshot().edges))))
    ("edges",
     serial_generate(Echo<double, CE>(sim.getHydroSnapshot().edges)))
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

/*! \brief Read dataset from file
  \param file File or group name
  \param field Name of field
  \param datatype Data type
  \return Dataset
*/
template<class FC> FC read_from_hdf5
(const H5::Group& file,
 const string& field,
 const H5::DataType& datatype)
{
  H5::DataSet dataset = file.openDataSet(field);
  H5::DataSpace filespace = dataset.getSpace();
  hsize_t dimes_out[2];
  filespace.getSimpleExtentDims(dimes_out);
  const size_t NX = static_cast<size_t>(dimes_out[0]);
  FC res;
  resize_if_necessary(res, NX);
  assert(NX == res.size());
  dataset.read(&res.front(), datatype);
  return res;
}

/*! \brief Combined lists of hydrodynamics variable to form computational cells
  \param density List of densities
  \param pressure List of pressures
  \param celerity List of proper velocities
  \return List of cells
 */
template<template<class> class CP> CP<Primitive> combine2cells
(const CP<double>& density,
 const CP<double>& pressure,
 const CP<double>& celerity)
{
  assert(density.size()==pressure.size() &&
	 density.size()==celerity.size());
  CP<Primitive> res;
  resize_if_necessary(res, density.size());
  //  array<CP<double>*, 3> master = {&density, &pressure, &celerity};
  for(size_t i=0;i<density.size();++i){
    res[i].Density = density[i];
    res[i].Pressure = pressure[i];
    res[i].Celerity = celerity[i];
    //    for(size_t j=0;j<3;++j)
    //      res[i][j] = master[j][i];
  }
  return res;
}

/*! \brief Read data from snapshot file
  \param fname File name
  \return Hydrodynamic snapshot
 */
template<template<class> class CE, template<class> class CP>
NewHydroSnapshot<CE, CP> read_hdf5_snapshot_new
(const string& fname)
{
  H5::H5File file(fname, H5F_ACC_RDONLY);
  const CE<double> edges = read_from_hdf5<CE<double> >
    (file, "edges", H5::PredType::NATIVE_DOUBLE);
  const CP<double> density = read_from_hdf5<CP<double> >
    (file, "density", H5::PredType::NATIVE_DOUBLE);
  const CP<double> pressure = read_from_hdf5<CP<double> >
    (file, "pressure", H5::PredType::NATIVE_DOUBLE);
  const CP<double> celerity = read_from_hdf5<CP<double> >
    (file, "celerity", H5::PredType::NATIVE_DOUBLE);
  assert(density.size()+1==edges.size());
  assert(density.size()==pressure.size());
  assert(density.size()==celerity.size());
  return NewHydroSnapshot<CE, CP>
    (edges,
     combine2cells<CP>
     (density,
      pressure,
      celerity));
}

#endif // HDF5_SNAPSHOT_HPP

