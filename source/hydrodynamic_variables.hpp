/*! \file hydrodynamic_variables.hpp
  \author Almog Yalinewich
  \brief Hydrodynamic variables
*/

#ifndef HYDRODYNAMIC_VARIABLES_HPP
#define HYDRODYNAMIC_VARIABLES_HPP 1

#include <string>
#include <vector>
#include <array>

using std::vector;
using std::string;
using std::array;

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

  Conserved(const array<double, 3>& source);

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

  NewConserved(const array<double, 3>& source);

  NewConserved(const NewConserved& source);

  NewConserved& operator=(const NewConserved& source);
};

//! \brief New primitive variables
class Primitive
{
public:

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

  //! \brief Density
  double Density;

  //! \brief Pressure
  double Pressure;

  //! \brief Celerity
  double Celerity;
};

/*! \brief Addition operator
  \param p1 Left argument
  \param p2 Right argument
  \return sum
 */
Primitive operator+(Primitive const& p1,
		    Primitive const& p2);

/*! \brief Subtraction operator
  \param p1 Left argument
  \param p2 Right argument
  \return Difference
 */
Primitive operator-(Primitive const& p1,
		    Primitive const& p2);

/*! \brief Division by a scalar
  \param p Primitive variable
  \param d Scalar
  \return ratio
 */
Primitive operator/(Primitive const& p,
		    double d);

/*! \brief Multiplication by a scalar
  \param d Scalar
  \param p Primitive variable
  \return Product
 */
Primitive operator*(double d,
		    Primitive const& p);

//! \brief Union for grid and primitive variables
class HydroSnapshot
{
public:

  /*! \brief Class constructor
    \param edges_i Grid
    \param cells_i List of primitive variables
   */
  HydroSnapshot(vector<double> const& edges_i,
		vector<Primitive> const& cells_i);

  /*! \brief Copy constructor
    \param source Source
   */
  HydroSnapshot(HydroSnapshot const& source);

  /*! \brief Copy constructor
    \param source Source
    return reference to self
  */
  HydroSnapshot& operator=(const HydroSnapshot& source);

  //! \brief Grid
  vector<double> edges;

  //! \brief Primitive variables
  vector<Primitive> cells;
};

/*! \brief Addition operator
  \param hs1 Left argument
  \param hs2 Right argument
  \return Sum
 */
HydroSnapshot operator+(HydroSnapshot const& hs1,
			HydroSnapshot const& hs2);

/*! \brief Subtraction operation
  \param hs1 Left argument
  \param hs2 Right argument
  \return Difference
 */
HydroSnapshot operator-(HydroSnapshot const& hs1,
			HydroSnapshot const& hs2);

/*! \brief Multiplication by a scalar
  \param d Scalar
  \param hs Hydrodynamic snapshot
  \return Product
 */
HydroSnapshot operator*(double d,
			HydroSnapshot const& hs);

/*! \brief Converts primitive variables to conserved (normalised by volume)
  \param p Primitive variables
  \return Conserved variables
*/
Conserved Primitive2Conserved_pv(Primitive const& p);

/*! \brief Addition of two vectors
  \param v1 First vector
  \param v2 Second vector
  \return Vector
*/
Conserved operator+(Conserved const& v1, 
		    Conserved const& v2);

/*! \brief Addition of two new conserved variables
  \param v1 First set of conserved variables
  \param v2 Second set of conserved variables
  \return Vector
*/
NewConserved operator+(NewConserved const& v1, 
		       NewConserved const& v2);

/*! \brief Addition of two vectors
  \param v1 First vector
  \param v2 Second vector
  \return Vector
*/
Conserved operator-(Conserved const& v1, Conserved const& v2);

/*! \brief Subtraction of two new conserved variables
  \param v1 First vector
  \param v2 Second vector
  \return Vector
*/
NewConserved operator-(NewConserved const& v1, NewConserved const& v2);

/*! \brief Subtraction of two vectors
  \param s Scalar
  \param v vector
  \return vector
*/
Conserved operator*(double s, Conserved const& v);

/*! \brief Subtraction of two vectors
  \param s Scalar
  \param v vector
  \return vector
*/
NewConserved operator*(double s, NewConserved const& v);

/*! \brief Scalar division
  \param s Scalar
  \param v vector
  \return vector
*/
Conserved operator/(Conserved const& v, double s);

#endif // HYDRODYNAMIC_VARIABLES_HPP
