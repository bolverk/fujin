/*! \file hydrodynamic_variables.hpp
  \author Almog Yalinewich
  \brief Hydrodynamic variables
*/

#ifndef HYDRODYNAMIC_VARIABLES_HPP
#define HYDRODYNAMIC_VARIABLES_HPP 1

#include <string>
#include <array>
#include <functional>

using std::string;
using std::array;
using std::function;
using std::pair;

//! \brief Conserved variables
class Conserved: public array<double, 3>
{
public:

  //! \brief Default constructor
  Conserved(void);

  /*! \brief Class constuctor
    \param mass Mass
    \param momentum Momentum
    \param energy Energy
  */
  Conserved(double mass,
	    double momentum,
	    double energy);

  explicit Conserved(const array<double, 3>& source);

  Conserved(const Conserved& source);

  Conserved& operator=(const Conserved& source);
  
  //! \brief Mass densty
  double& Mass;

  //! \brief Momentum density
  double& Momentum;

  //! \brief Energy density
  double& Energy;
};

//! \brief New conserved variable. These will not become degenerate at ultra relativistic flows
class NewConserved: public array<double, 3>
{
public:

  //! \brief Mass density
  double& mass;

  //! \brief Sum of energy and momentum densities
  double& positive;

  //! \brief Difference of energy and momentum densities
  double& negative;

  //! \brief Null constructor. Sets everything to zero.
  NewConserved(void);

  /*! \brief Class constructor
    \param mass_i Mass density
    \param positive_i Sum of energy and momentum densities
    \param negative_i Difference of energy and momentum densities
   */
  NewConserved(double mass_i,
	       double positive_i,
	       double negative_i);

  explicit NewConserved(const array<double, 3>& source);

  NewConserved(const NewConserved& source);

  NewConserved& operator=(const NewConserved& source);
};

//! \brief New primitive variables
class Primitive: public array<double, 3>
{
public:

  explicit Primitive(const array<double, 3>& source);

  //! \brief Null constructor
  Primitive(void);

  /*! \brief Class constructor
    \param density Density
    \param pressure Pressure
    \param celerity Celerity
  */
  Primitive(double density,
	    double pressure,
	    double celerity);

  Primitive(const Primitive& source);

  Primitive& operator=(const Primitive& source);

  //! \brief Density
  double& Density;

  //! \brief Pressure
  double& Pressure;

  //! \brief Celerity
  double& Celerity;
};

/*! \brief Converts primitive variables to conserved (normalised by volume)
  \param p Primitive variables
  \return Conserved variables
*/
Conserved Primitive2Conserved_pv(Primitive const& p);

template<class T> T bin_op
(const T& p1,
 const T& p2,
 function<double(double, double)> func)
{
  T res;
  transform(p1.begin(),
	    p1.end(),
	    p2.begin(),
	    res.begin(),
	    func);
  return res;
}

template<class T> T une_op
(const T& t,
 function<double(double)> func)
{
  T res;
  transform(t.begin(),
	    t.end(),
	    res.begin(),
	    func);
  return res;
}

template<class T> typename std::enable_if<std::is_base_of<array<double,3>, T>::value, T>::type operator+
(const T& t1,
 const T& t2)
{
  return bin_op(t1, t2, std::plus<double>());
}

template<class T> typename std::enable_if<std::is_base_of<array<double,3>, T>::value, T>::type operator-
(const T& t1,
 const T& t2)
{
  return bin_op(t1, t2, std::minus<double>());
}

template<class T> typename std::enable_if<std::is_base_of<array<double, 3>, T>::value, T>::type operator*
(const T& t,
 double s)
{
  return une_op(t, [&s](double d){return s*d;});
}

template<class T> typename std::enable_if<std::is_base_of<array<double, 3>, T>::value, T>::type operator*
(double s,
 const T& t)
{
  return une_op(t, [&s](double d){return s*d;});
}

template<class T> typename std::enable_if<std::is_base_of<array<double, 3>, T>::value, T>::type operator/
(const T& t,
 double s)
{
  return une_op(t, [&s](double d){return d/s;});
}

template<template<class> class CE, template<class> class CP> class NewHydroSnapshot: public pair<CE<double>, CP<Primitive> >
{
public:

  NewHydroSnapshot(const CE<double>& edges_i,
		   const CP<Primitive>& cells_i):
    pair<CE<double>, CP<Primitive> >(edges_i, cells_i),
    edges((*this).first),
    cells((*this).second) {}

  NewHydroSnapshot(const NewHydroSnapshot<CE, CP>& source):
    NewHydroSnapshot(source.edges, source.cells) {}

  //  NewHydroSnapshot(const HydroSnapshot& source):
  //    NewHydroSnapshot(source.edges, source.cells) {}

  NewHydroSnapshot& operator=(const NewHydroSnapshot& source)
  {
    pair<CE<double>, CP<Primitive> >::operator=(source);
    return *this;
  }

  CE<double>& edges;
  CP<Primitive>& cells;  
};

#endif // HYDRODYNAMIC_VARIABLES_HPP
