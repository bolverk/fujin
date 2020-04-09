#include "hdf5_snapshot.hpp"
#include "hdf5_utils.hpp"
#include "diagnostics.hpp"
#include <assert.h>

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
