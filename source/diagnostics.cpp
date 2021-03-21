#include <fstream>
#include <boost/foreach.hpp>
#include "diagnostics.hpp"
#include "advanced_hydrodynamic_variables.hpp"
#include "vector_initializer.hpp"
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
  /*! \brief Check if all members are true
    \param flags Lazy list of booleans
    \return true if all are true, false otherwise
   */
  bool all_true(const Index2Member<bool>& flags)
  {
    for(size_t i=0;i<flags.getLength();++i){
      if(!flags(i))
	return false;
    }
    return true;
  }
}

namespace {

  //! \brief Checks if conserved variables and primitives are thermodynamically consistent
#if SCAFFOLDING != 1
  template<class CE, class CP>
#endif // SCAFFOLDING
  class ConservedPrimitiveConsistencyChecker: public Index2Member<bool>
  {
  public:

    /*! \brief Class constructor
      \param sim Hydrodynamic simulation
      \param thres Threshold
     */
    ConservedPrimitiveConsistencyChecker
    (
#if SCAFFOLDING == 1
     const SRHDSimulation& sim,
#else
     const SRHDSimulation<CE, CP>& sim,
#endif // SCAFFOLDING
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
#if SCAFFOLDING == 1
    const SRHDSimulation& sim_;
#else
    const SRHDSimulation<CE, CP>& sim_;
#endif // SCAFFOLDING
    //! \brief Threshold
    const double thres_;
  };
}

bool ConservedPrimitiveConsistency
(SRHDSimulation const& sim, 
 double thres)
{
  return all_true(ConservedPrimitiveConsistencyChecker(sim,thres));
}

namespace {

  //! \brief Calculates the stresses
  class StressCalculator: public Index2Member<double>
  {
  public:

    /*! \brief Class constructor
      \param sim Simulation
      \param idx Index to member
     */
    StressCalculator(const SRHDSimulation& sim,
		     size_t idx):
      sim_(sim), idx_(idx) {}

    size_t getLength(void) const
    {
      return sim_.getHydroSnapshot().cells.size();
    }

    double operator()(size_t i) const
    {
      return sim_.getRestMasses().at(i)*sim_.getConserved().at(i)[idx_];
    }

  private:
    //! \brief Simulation
    SRHDSimulation const& sim_;
    //! \brief Pointer to member
    size_t idx_;
    //    double NewConserved::* pcm_;
  };

  class CellVolumes: public Index2Member<double>
  {
  public:

    explicit CellVolumes(const SRHDSimulation& sim):
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
    const SRHDSimulation& sim_;
  };

  template<class T> class ElementwiseSum:
    public Index2Member<T>
  {
  public:

    ElementwiseSum
    (const Index2Member<T>& list_1,
     const Index2Member<T>& list_2):
      list_1_(list_1),
      list_2_(list_2)
    {
      assert(list_1.getLength()==list_2.getLength());
    }

    size_t getLength(void) const
    {
      return list_1_.getLength();
    }

    T operator()(size_t i) const
    {
      return list_1_(i) + list_2_(i);
    }

  private:
    const Index2Member<T>& list_1_;
    const Index2Member<T>& list_2_;
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

double TotalEnergy(SRHDSimulation const& sim)
{
  return 0.5*
    (sum_all
     (ElementwiseProduct<double,double,double>
      (CellVolumes(sim),
       ElementwiseSum<double>
       (StressCalculator(sim,1),
	StressCalculator(sim,2)))));
  //  return 0.5*(sum_all(StressCalculator(sim,&NewConserved::positive))+
  //	      sum_all(StressCalculator(sim,&NewConserved::negative)));
}

double TotalMomentum(SRHDSimulation const& sim)
{
  return 0.5*(sum_all(StressCalculator(sim,1))-
	      sum_all(StressCalculator(sim,2)));
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

bool VerticesIncreasingOrder(SRHDSimulation const& sim)
{
  return all_true
    (IsIncreasing<double>(Echo<double>(sim.getHydroSnapshot().edges)));
}

namespace {

  //! \brief Retrieves the centres of cells
  #if 0
  class CellCenterGetter: public Index2Member<double>
  {
  public:

    /*! \brief Class constructor
      \param sim Simulation
     */
    explicit CellCenterGetter(const SRHDSimulation& sim):
      edges_(sim.getHydroSnapshot().edges) {}

    size_t getLength(void) const
    {
      return edges_.size()-1;
    }

    double operator()(size_t i) const
    {
      return 0.5*(edges_[i]+edges_[i+1]);
    }

  private:
    //! \brief Grid
    const vector<double> edges_;
  };
#endif // 0
}

PrimitivePropertyGetter::PrimitivePropertyGetter
(const SRHDSimulation& sim,
 size_t idx):
  cells_(sim.getHydroSnapshot().cells), idx_(idx) {}

size_t PrimitivePropertyGetter::getLength(void) const
{
  return cells_.size();
}

double PrimitivePropertyGetter::operator()(size_t i) const
{
  return cells_[i][idx_];
}

/*
void write_snapshot(SRHDSimulation const& sim,
		    string const& fname,
		    int precision)
{
  vector<unique_ptr<Index2Member<double> > > properties;
  properties.push_back(make_unique<CellCenterGetter>(sim));
  for(size_t i=0;i<3;++i)
    properties.push_back(make_unique<PrimitivePropertyGetter>(sim, i));
  std::ofstream f(fname.c_str());
  f.precision(precision);
  for(size_t i=0;i<sim.getHydroSnapshot().cells.size();++i){
    for(auto &p : properties)
      f << (*p)(i) << " ";
    f << "\n";
  }
  f.close();
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
