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

  /*! \brief Copy constructor
    \param source Source
   */
  explicit Conserved(const array<double, 3>& source);

  /*! \brief Copy constructor
    \param source Source
   */
  Conserved(const Conserved& source);

  /*! \brief Copy assignment
    \param source Source
    \return Reference to self
   */
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

  /*! \brief Copy constructor
    \param source Source
   */
  explicit NewConserved(const array<double, 3>& source);

  /*! \brief Copy constructor
    \param source Source
   */
  NewConserved(const NewConserved& source);

  /*! \brief Assignment operator
    \param source Source
    \return Reference to self
   */
  NewConserved& operator=(const NewConserved& source);
};

//! \brief New primitive variables
class Primitive: public array<double, 3>
{
public:

  /*! \brief Class constructor
    \param source Values of primitive variables
   */
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

  /*! \brief Copy constructor
    \param source Source
   */
  Primitive(const Primitive& source);

  /*! \brief Assignment operator
    \param source Source
    \return Reference to self
   */ 
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

/*! \brief Binary operation
  \param p1 First argument
  \param p2 Second argument
  \param func Binary function
  \return Result of func applied to all member
 */
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

/*! \brief Unary operation
  \param t Argument
  \param func Unary functions
  \return Result of unary function
 */
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

/*! \brief Subtraction operator
  \param t1 First argument
  \param t2 Second argument
  \return The result of t2 subtracted from t1
 */
template<class T> typename std::enable_if<std::is_base_of<array<double,3>, T>::value, T>::type operator+
(const T& t1,
 const T& t2)
{
  return bin_op(t1, t2, std::plus<double>());
}

/*! \brief Subtraction operator
  \param t1 First argument
  \param t2 Second argument
  \return The difference of t2 from t1
 */
template<class T> typename std::enable_if<std::is_base_of<array<double,3>, T>::value, T>::type operator-
(const T& t1,
 const T& t2)
{
  return bin_op(t1, t2, std::minus<double>());
}

/*! \brief Multiplication by a scalar
  \param t Argument
  \param s Scalar
  \return The product t*s
 */
template<class T> typename std::enable_if<std::is_base_of<array<double, 3>, T>::value, T>::type operator*
(const T& t,
 double s)
{
  return une_op(t, [&s](double d){return s*d;});
}

/*! \brief Multiplication by scalar
  \param s Scalar
  \param t Argument
  \return The product t*s
 */
template<class T> typename std::enable_if<std::is_base_of<array<double, 3>, T>::value, T>::type operator*
(double s,
 const T& t)
{
  return une_op(t, [&s](double d){return s*d;});
}

/*! \brief Division by a scalar
  \param t Argument
  \param s Scalar
  \return The product (1/s)*t;
 */
template<class T> typename std::enable_if<std::is_base_of<array<double, 3>, T>::value, T>::type operator/
(const T& t,
 double s)
{
  return une_op(t, [&s](double d){return d/s;});
}

//! \brief Hydrodynamic snapshot
template<template<class> class CE, template<class> class CP> class NewHydroSnapshot: public pair<CE<double>, CP<Primitive> >
{
public:

  /*! \brief Class constructor
    \param edges_i Cell edges
    \param cells_i Computational cells
   */
  NewHydroSnapshot(const CE<double>& edges_i,
		   const CP<Primitive>& cells_i):
    pair<CE<double>, CP<Primitive> >(edges_i, cells_i),
    edges((*this).first),
    cells((*this).second) {}

  /*! \brief Copy constructor
    \param source Source
   */
  NewHydroSnapshot(const NewHydroSnapshot<CE, CP>& source):
    NewHydroSnapshot(source.edges, source.cells) {}

  /*! \brief Assignment operator
    \param source Source
    \return Reference to self
   */
  NewHydroSnapshot& operator=(const NewHydroSnapshot& source)
  {
    pair<CE<double>, CP<Primitive> >::operator=(source);
    return *this;
  }

  //! \brief Cell edges
  CE<double>& edges;
  //! \brief Computational cells
  CP<Primitive>& cells;  
};

#endif // HYDRODYNAMIC_VARIABLES_HPP
