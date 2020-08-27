#include "hdf5_snapshot.hpp"
#include "hdf5_utils.hpp"
#include "diagnostics.hpp"
#include <assert.h>
#include <H5Cpp.h>
#include <map>
#include <boost/foreach.hpp>

using namespace H5;

namespace {

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
}

void write_hdf5_snapshot(const SRHDSimulation& sim,
			 const string& fname)
{
  (HDF5Shortcut(fname))
    ("position",
     serial_generate(MidValues(Echo<double>(sim.getHydroSnapshot().edges))))
	("edges",
     serial_generate(Echo<double>(sim.getHydroSnapshot().edges)))
    ("density",
     serial_generate(PrimitivePropertyGetter(sim,&Primitive::Density)))
    ("pressure",
     serial_generate(PrimitivePropertyGetter(sim,&Primitive::Pressure)))
    ("celerity",
     serial_generate(PrimitivePropertyGetter(sim,&Primitive::Celerity)))
    ("time",
     vector<double>(1,sim.GetTime()))
    ("cycle",
     vector<double>(1,sim.GetCycle()));
}

namespace {
  template<class T> vector<T> read_vector_from_hdf5
  (const Group& file,
   const string& field,
   const DataType& datatype)
  {
    DataSet dataset = file.openDataSet(field);
    DataSpace filespace = dataset.getSpace();
    hsize_t dimes_out[2];
    filespace.getSimpleExtentDims(dimes_out);
    const size_t NX = static_cast<size_t>(dimes_out[0]);
    vector<T> result(NX);
    dataset.read(&result[0], datatype);
    return result;
  }

  vector<double> read_double_vector_from_hdf5
  (const Group& file,
   const string& field)
  {
    return read_vector_from_hdf5<double>
      (file, field, PredType::NATIVE_DOUBLE);
  }

  vector<Primitive> combine2cells
  (const vector<double>& density,
   const vector<double>& pressure,
   const vector<double>& celerity)
  {
    assert(density.size()==pressure.size() &&
	   density.size()==celerity.size());
    vector<Primitive> res(density.size());
    for(size_t i=0;i<density.size();++i)
      res[i] = Primitive(density[i],
			 pressure[i],
			 celerity[i]);
    return res;
  }
}

HydroSnapshot read_hdf5_snapshot
(const string& fname)
{
  H5File file(fname, H5F_ACC_RDONLY);
  const vector<string> field_list =
    {"edges",
     "density",
     "pressure",
     "celerity"};
  std::map<string,vector<double> > data;
  BOOST_FOREACH(string field, field_list)
    {
      data[field] = read_double_vector_from_hdf5(file, field);
    }
  assert(data["density"].size()+1==data["edges"].size());
  assert(data["density"].size()==data["pressure"].size());
  assert(data["density"].size()==data["celerity"].size());
  return HydroSnapshot(data["edges"],
		       combine2cells(data["density"],
				     data["pressure"],
				     data["celerity"]));
}
