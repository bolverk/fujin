/*! \file diagnostics.hpp
  \brief Diagnostic functions for srhd_sim
  \author Almog Yalinewich
*/

#ifndef DIAGNOSTICS_HPP
#define DIAGNOSTICS_HPP 1

#include <fstream>
#include "hydrodynamic_variables.hpp"
#include "srhd_sim.hpp"

/*! \brief Approximate comparison
  \details Verifies that two values are very close to one another
  \param v1 First value
  \param v2 Second value
  \param thres Threshold
  \return True if the difference between v1 and v2 is below thres, false otherwise
 */
bool ApproxCompare(double v1, double v2, double thres);

/*! \brief Verifies that the primitive variables are consistent with the conserved variable, up to a specified accuracy
  \param p Primitive variables
  \param c Conserved variables
  \param eos Equation of state
  \param thres Threshold
  \return True if both parameters agree, false otherwise
 */
bool ConservedPrimitiveConsistency(Primitive const& p,
				   NewConserved const& c,
				   const EquationOfState& eos,
				   double thres);

/*! \brief Verifies that the primitive variables in all the cells are consistent with the conserved variables
  \param sim Hydrodynamic simulation
  \param thres Threshold
  \return True if all cells are consistent, false otherwise
 */
template<template<class> class CE, template<class> class CP>
bool ConservedPrimitiveConsistency
(const SRHDSimulation<CE, CP>& sim, 
 double thres);

template<template<class> class CE, template<class> class CP> class CellVolumes: public Index2Member<double>
{
public:

  explicit CellVolumes(const SRHDSimulation<CE, CP>& sim):
    sim_(sim) {}

  size_t getLength(void) const override
  {
    return sim_.getHydroSnapshot().cells.size();
  }

  double operator()(size_t i) const override
  {
    return sim_.getVolume(i);
  }

private:

  const SRHDSimulation<CE, CP>& sim_;
};

template<template<class> class CE, template<class> class CP> class StressCalculator: public Index2Member<double>
{
public:

  StressCalculator(const SRHDSimulation<CE, CP>& sim,
		   size_t idx):
    sim_(sim), idx_(idx) {}

  size_t getLength(void) const override
  {
    return sim_.getHydroSnapshot().cells.size();
  }

  double operator()(size_t i) const override
  {
    return sim_.getRestMasses().at(i)*sim_.getConserved().at(i)[idx_];
  }

private:

  const SRHDSimulation<CE, CP>& sim_;
  const size_t idx_;
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

    size_t getLength(void) const override
    {
      return list_1_.getLength();
    }

    T operator()(size_t i) const override
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

    size_t getLength(void) const override
    {
      return list_1_.getLength();
    }

    T3 operator()(size_t i) const override
    {
      return list_1_(i)*list_2_(i);
    }

  private:
    const Index2Member<T1>& list_1_;
    const Index2Member<T2>& list_2_;
  };

/*! \brief Returns the total energy
  \param sim Hydrodynamic simulation
  \return Total energy
 */
template<template<class> class CE, template<class> class CP>
double TotalEnergy(const SRHDSimulation<CE, CP>& sim)
{
    return 0.5*
    (sum_all
     (ElementwiseProduct<double,double,double>
      (
      CellVolumes<CE,CP>(sim),
       ElementwiseSum<double>
      (StressCalculator<CE, CP>(sim,1),
       StressCalculator<CE, CP>(sim,2)))));
}

/*! \brief Returns the total momentum
  \param sim Hydrodynamic simulation
  \return Total momentum
 */
template<template<class> class CE, template<class> class CP>
double TotalMomentum(SRHDSimulation<CE, CP> const& sim)
{
  return 0.5*(sum_all(StressCalculator<CE, CP>(sim,1))-
	      sum_all(StressCalculator<CE, CP>(sim,2)));
}

/*! \brief Checks that the vertices are in increasing order
  \param sim Reference to simulation
  \return True if they are in increasing order, false otherwise
 */
template<template<class> class CE, template<class> class CP>
bool VerticesIncreasingOrder(const SRHDSimulation<CE, CP>& sim);

/*! \brief Writes list of data to file as ascii
  \param v Data
  \param fname Name of file
  \param precision Number of digits
 */
template <class T> void write_to_file(vector<T> const& v, string const& fname,
				      int precision = 6)
{
  std::ofstream f(fname.c_str());
  f.precision(precision);
  for_each(v.begin(),
	   v.end(),
	   [&f](const T& t)
	   {f << t << "\n";});
  f.close();
}

/*! \brief Writes grid and cells to a file
  \param sim Hydrodynamic simulation
  \param fname Name of file
  \param precision Number of digits
 */
/*
void write_snapshot(SRHDSimulation const& sim,
		    string const& fname,
		    int precision=6);
*/

/*! \brief Writes a single number to a file
  \param num Number
  \param fname Name of file
  \param precision Number of digits to write
 */
void write_number(double num,
		  string const& fname,
		  int precision=14);

//! \brief Retrieves properties from a list of primitive variables
template<template<class> class CE, template<class> class CP>
class PrimitivePropertyGetter: public Index2Member<double>
{
public:

  /*! \brief Class constructor
    \param sim Simulation
    \param idx Index to member
  */
  PrimitivePropertyGetter
  (const SRHDSimulation<CE, CP>& sim,
   size_t idx):
  cells_(sim.getHydroSnapshot().cells),
    idx_(idx) {}    

  size_t getLength(void) const override
  {
    return cells_.size();
  }

  double operator()(size_t i) const override
  {
    return cells_[i][idx_];
  }

private:
  //! \brief List of primitive variables
  const CP<Primitive>& cells_;

  //! \brief Pointer to member of primitives
  const size_t idx_;
};

/*! \brief Writes all vector member to a file
  \param list Data to be written
  \param fname Name of file
  \param precision Precision
 */
/*
void write_vector(const vector<double>& list,
		  const string& fname,
		  int precision=14);
*/

#endif // DIAGNOSTICS_HPP
