#include <fstream>
#include <boost/foreach.hpp>
#include "diagnostics.hpp"
#include "advanced_hydrodynamic_variables.hpp"
#include <assert.h>
#include <memory>

using std::unique_ptr;
using std::make_unique;

bool ApproxCompare(double v1, double v2, double thres)
{
  const double eps = 1e-12;
  return (fabs(v1-v2)/(eps+fabs(v1)+fabs(v2)))<thres;
}

bool ConservedPrimitiveConsistency(Primitive const& p,
				   NewConserved const& c,
				   const EquationOfState& eos,
				   double thres)
{
  NewConserved ref = primitive_to_new_conserved(p, eos);
  return ApproxCompare(ref.mass,c.mass,thres)&&
    ApproxCompare(ref.positive,c.positive,thres)&&
    ApproxCompare(ref.negative,c.negative,thres);
}

namespace {
  /*
  static bool all_true(const Index2Member<bool>& flags)
  {
    for(size_t i=0;i<flags.getLength();++i){
      if(!flags(i))
	return false;
    }
    return true;
  }
  */
}

namespace {

  //! \brief Checks if conserved variables and primitives are thermodynamically consistent
  template<template<class> class CE, template<class> class CP>
  class ConservedPrimitiveConsistencyChecker: public Index2Member<bool>
  {
  public:

    /*! \brief Class constructor
      \param sim Hydrodynamic simulation
      \param thres Threshold
    */
    ConservedPrimitiveConsistencyChecker
    (
     const SRHDSimulation<CE, CP>& sim,
     double thres):
      sim_(sim), thres_(thres) {}

    size_t getLength(void) const
    {
      return sim_.getHydroSnapshot().cells.size();
    }

    bool operator()(size_t i) const
    {
      return ConservedPrimitiveConsistency
	(sim_.getHydroSnapshot().cells[i],
	 sim_.getConserved().at(i),
	 sim_.getEOS(),
	 thres_);
    }

  private:
    //! \brief Simulation
    const SRHDSimulation<CE, CP>& sim_;
    //! \brief Threshold
    const double thres_;
  };
}

/*
template<class CE, class CP>
bool ConservedPrimitiveConsistency
(
 const SRHDSimulation<CE, CP>& sim,
 double thres)
{
  return all_true(ConservedPrimitiveConsistencyChecker<CE,CP>(sim,thres));
}
*/

namespace {

  template<template<class> class CE, template<class> class CP>
  class CellVolumes: public Index2Member<double>
  {
  public:

    explicit CellVolumes
    (const SRHDSimulation<CE, CP>& sim):
      sim_(sim) {}

    size_t getLength(void) const
    {
      return sim_.getHydroSnapshot().cells.size();
    }

    double operator()(size_t i) const
    {
      return sim_.getVolume(i);	
    }

  private:
    const SRHDSimulation<CE, CP>& sim_;
  };

  template<class T1, class T2, class T3>
  class ElementwiseProduct:
    public Index2Member<T3>
  {
  public:

    ElementwiseProduct
    (const Index2Member<T1>& list_1,
     const Index2Member<T2>& list_2):
      list_1_(list_1),
      list_2_(list_2)
    {
      assert(list_1.getLength()==list_2.getLength());
    }

    size_t getLength(void) const
    {
      return list_1_.getLength();
    }

    T3 operator()(size_t i) const
    {
      return list_1_(i)*list_2_(i);
    }

  private:
    const Index2Member<T1>& list_1_;
    const Index2Member<T2>& list_2_;
  };
}

//! \brief Auxiliary class for checking whether a list is increasing
template<class T> class IsIncreasing: public Index2Member<bool>
{
public:

  /*! \brief Class constructor
    \param i2m Lazy list
  */
  explicit IsIncreasing(const Index2Member<T>& i2m):
    i2m_(i2m) {}

  size_t getLength(void) const
  {
    return i2m_.getLength()-1;
  }

  bool operator()(size_t i) const
  {
    return i2m_(i+1)>i2m_(i);
  }

private:

  //! \brief Lazy list
  const Index2Member<T>& i2m_;
};

/*
template<class CE, class CP>
bool VerticesIncreasingOrder(const SRHDSimulation<CE, CP>& sim)
{
  return all_true
    (IsIncreasing<double>(Echo<double>(sim.getHydroSnapshot().edges)));
}
*/

void write_number(double num,
		  string const& fname,
		  int precision)
{
  std::ofstream f(fname.c_str());
  f.precision(precision);
  f << num << "\n";
  f.close();
}

/*
  void write_vector(const vector<double>& list,
  const string& fname,
  int precision)
  {
  std::ofstream f(fname.c_str());
  f.precision(precision);
  BOOST_FOREACH(double itm, list)
  {
  f << itm << "\n";
  }
  f.close();
  }
*/
