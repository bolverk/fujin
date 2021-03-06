#include "hdf5_utils.hpp"
#include <H5Cpp.h>
#include <cassert>
#include <algorithm>
//#include "utils.hpp"

using namespace H5;
using std::transform;

HDF5Shortcut::HDF5Shortcut(const string& fname):
  fname_(fname), double_data_(), int_data_() {}

HDF5Shortcut& HDF5Shortcut::operator()(const string& field_name,
				       const vector<double>& data_set)
{
  assert(!data_set.empty());
  double_data_.push_back(pair<string,vector<double> >
			 (field_name,data_set));
  return *this;
}

HDF5Shortcut& HDF5Shortcut::operator()(const string& field_name,
				       const vector<int>& data_set)
{
  int_data_.push_back(pair<string,vector<int> >
		      (field_name,data_set));
  return *this;
}

namespace {

  /*! \brief Performs static cast for all items in a list
    \param source Source list
    \return source list cast to new type
   */
  template<class T, class S> vector<T> list_static_cast(const vector<S>& source)
  {
    vector<T> res(source.size());
    transform(source.begin(),
	      source.end(),
	      res.begin(),
	      [](const S& s)
	      {return static_cast<T>(s);});
    return res;
  }

  /*! \brief Writes a list of integers as an HDF5 dataset
    \param file HDF5 file
    \param num_list List of integers
    \param caption Name of dataset
   */
  void write_std_vector_to_hdf5
  (H5File& file,
   vector<int> const& num_list,
   string const& caption)
  {
    hsize_t dimsf[1];
    dimsf[0] = static_cast<hsize_t>(num_list.size());
    DataSpace dataspace(1, dimsf);

    IntType datatype(PredType::NATIVE_UINT);
    datatype.setOrder(H5T_ORDER_LE);

    // Modify dataset creation property to enable chunking
    DSetCreatPropList  plist;
    if(dimsf[0]>100000)
      dimsf[0]=100000;
    plist.setChunk(1,dimsf);
    plist.setDeflate(6);

    DataSet dataset = file.createDataSet(H5std_string(caption),
					 datatype,
					 dataspace,plist);

    const vector<unsigned> buf = list_static_cast<unsigned,int>(num_list);
    dataset.write(&buf[0],PredType::NATIVE_UINT);
  }

  /*! \brief Writes a vector to an hdf5 dataset
    \param file HDF5 file
    \param num_list List of numbers
    \param caption Name of dataset
   */
  void write_std_vector_to_hdf5
  (H5File& file,
   vector<double> const& num_list,
   string const& caption)
  {
    hsize_t dimsf[1];
    dimsf[0] = static_cast<hsize_t>(num_list.size());
    DataSpace dataspace(1, dimsf);

    FloatType datatype(PredType::NATIVE_DOUBLE);
    datatype.setOrder(H5T_ORDER_LE);

    // Modify dataset creation property to enable chunking
    DSetCreatPropList  plist;
    if(dimsf[0]>100000)
      dimsf[0]=100000;
    plist.setChunk(1,dimsf);
    plist.setDeflate(6);

    DataSet dataset = file.createDataSet(H5std_string(caption),
					 datatype,
					 dataspace,plist);

    dataset.write(&num_list[0],PredType::NATIVE_DOUBLE);
  }
}

HDF5Shortcut::~HDF5Shortcut(void)
{
  H5File file(H5std_string(fname_), H5F_ACC_TRUNC);
  for(size_t i=0;i<double_data_.size();++i)
    write_std_vector_to_hdf5(file, double_data_[i].second, double_data_[i].first);
  for(size_t i=0;i<int_data_.size();++i)
    write_std_vector_to_hdf5(file, int_data_[i].second, int_data_[i].first);
}
