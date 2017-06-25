//! \file igrs.hpp
//! \author Almog Yalinewich
//! \brief Riemann Solver for a ideal reletavistic photon gas

#ifndef IGRS_HPP
#define IGRS_HPP

#include "utilities.hpp"

/*! \brief Calculates the dimensionless velocity along the Hugoniot curve
  \param p2 {New pressure} \param p1 Old pressure
  \param g Adiabatic index
  \return Dimensionless velocity
*/
double HugoniotVelocity(double p2, double p1, double g);

/*! \brief Derivative of HugoniotVelocity with respect to p2
	\param p2 New pressure
	\param p1 Old pressure
	\param g Adiabatic inde
	\return Derivative of the dimensionless velocity with respect to pressure
*/
double dHugoniotVelocity(double p2, double p1, double g);

/*! \brief Calculates the dimensionless velocity along the isentrope curve
  \param p2 New pressure
  \param p1 Old pressure
  \param g Adiabatc index
  \return Dimensionless velocity
 */
double IsentropeVelocity(double p2, double p1, double g);

/*! \brief Calculates the derivative of IsentropeVelocity with respect to p2
  \param p2 New pressure
  \param p1 Old pressure
  \param g Adiabatic index
  \return Derivative of the dimensionless velocity with respect to pressure
 */
double dIsentropeVelocity(double p2, double p1, double g);

/*! \brief Calculates the dimensionless velocity along the hydrodynamic curve
  \details If p2>p1, then this curve is the hugoniot, otherwise, it is the isentrope
  \param p2 New pressure
  \param p1 Old pressure
  \param g Adiabatic index
  \return Dimensionless velocity
 */
double HydrodynamicVelocity(double p2, double p1, double g);

/*! \brief Calculates the derivative of HydrodynamicVelocity with respect to p2
  \param p2 New pressure
  \param p1 Old pressure
  \param g Adiabatic index
  \return Derivative of the dimensionless velocity with respect to pressure
 */
double dHydrodynamicVelocity(double p2, double p1, double g);

/*! \brief Acoustic approximation to the pressure at the interface of the Riemann problem
  \param pl Left pressure
  \param pr Right pressure
  \param bl Left dimensionless velocity
  \param br Right dimensionless velocity
  \param g Adiabatic index
  \return Pressure
 */
 
double PsAcoustic(double pl, double pr, double bl, double br, double g);
/*! \brief Calculates the dimensionless velocity of the left (negative) side in the Riemann problem
  \param b0 Initial dimensionless velocity
  \param p New pressure
  \param p0 Old pressure
  \param g Adiabatic index
  \return Dimensionless velocity
 */
double LeftVelocity(double b0, double p, double p0, double g);

/*! \brief Calculates the derivative of LeftVelocity with respect to p
  \param b0 Initial dimensionless velocity
  \param p New pressure
  \param p0 Old pressure
  \param g Adiabatic index
  \return Derivative of the dimensionless velocity with respect to pressure
 */
double dLeftVelocity(double b0, double p, double p0, double g);

/*! \brief Calculates the dimensionless velocity on the right (positive) side in the Riemann problem
  \param b0 Initial dimensionless velocity
  \param p New pressure
  \param p0 Old pressure
  \param g Adiabatic index
  \return Dimensionless velocity
 */
double RightVelocity(double b0, double p, double p0, double g);

/*! \brief Calculates the derivative of RightVelocity with respect to p
  \param b0 Initial dimensionless velocity
  \param p New pressure
  \param p0 Old pressure
  \param g Adiabatic index
  \return Derivative of the dimensionless velocity with respect to pressure
 */
double dRightVelocity(double b0, double p, double p0, double g);

/*! \brief Evaluates the transcendental equation of the Riemann problem
  \param bl Left dimensionless velocity
  \param br Right dimensionless velocity
  \param pl Left pressure
  \param pr Right pressure
  \param p New pressure
  \param g Adiabatic index
  \return Double
 */
double TranscendentalEquation(double bl, double br, double pl, 
			      double pr, double p, double g);
/*! \brief Calculates the derivative of TranscendentalEquation with respect to p
  \param bl Left dimensionless velocity
  \param br Right dimensionless velocity
  \param pl Left pressure
  \param pr Right pressure
  \param p New pressure
  \param g Adiabatic index
  \return Derivative of the equation
 */
double dTranscendentalEquation(double bl, double br, double pl, 
			       double pr, double p, double g);

/*! \brief Interface for the function TranscendentalEquation
 */
class iTranscendentalEquation: public SVDifferentiable
{
private:
  //! \brief Left dimensionless velocity
  const double  bl;
  //! \brief Right dimensionless velocity
  const double br;
  //! \brief Left pressure
  const double pl;
  //! \brief Right pressure
  const double  pr;
  //! \brief Adiabatic index
  const double g;
public:
  /*! \brief Class costructor
    \param bli Left dimensionless velocity
    \param bri Right dimensionless velocity
    \param pli Left pressure
    \param pri Right pressure
    \param gi Adiabatic index
   */
  iTranscendentalEquation(double bli, double bri,
			  double pli, double pri,
			  double gi);
  /*! \brief Evaluates the Transcendental equation
    \param x Pressure
    \return Double
  */
  double operator()(double x) const;

  /*! \brief Evaluates the derivative to the transcendental equation with respect to pressure
    \param x Pressure
    \return Double
   */
  double diff(double x) const;
};

/*! \brief Solution to the Riemann problem
 */
struct RiemannSolution
{
  //! \brief Pressure at the interface
  double Pressure;
  //! \brief Velocity at the interface
  double Velocity;
};

/*! \brief Solves the transcendental equation for the pressure
  \param bl Left dimensionless velocity
  \param br Right dimensionless velcoity
  \param pl Left pressure
  \param pr Right pressure
  \param g Adiabatic index
  \return Pressure
 */
double SolveTranscendentalEquation(double bl, double br, 
				   double pl, double pr, 
				   double g);

/*! \brief Solves the Riemann problem
  \details Return the pressure and velocity at the interface
  \param bl Left dimensionless velocity
  \param br Right dimensionless velcoity
  \param pl Left pressure
  \param pr Right pressure
  \param g Adiabatic index
  \return Pressure and velocity
 */
RiemannSolution RiemannSolve(double bl, double br,
			     double pl, double pr,
			     double g);

#endif
