//! \file igrs.cpp
//! \author Almog Yalinewich
//! \brief Riemann solver for an ideal, relativistic gas

#include <math.h>
#include "igrs.hpp"
#include "utilities.hpp"

double HugoniotVelocity(double p2, double p1, double g)
{
  return (-p1 + p2)/(sqrt(p1/(-1 + g) + p2)*sqrt((-1 + g)*p1 + p2));
}

double dHugoniotVelocity(double p2, double p1, double g)
{
  return pow(g,2)*p1*(p1+p2)/
    (2*sqrt(p2+p1/(g-1))*
     pow(p2+p1*(g-1),1.5)*
     (p1+(g-1)*p2));
}

double IsentropeVelocity(double p2, double p1, double g)
{
  return tanh(sqrt(g-1)*log(p2/p1)/g);
}

double dIsentropeVelocity(double p2, double p1, double g)
{
  return sqrt(g-1)*pow(cosh(sqrt(g-1)*log(p2/p1)/g),-2)/g/p2;
}

double HydrodynamicVelocity(double p2, double p1, double g)
{
  if (p2>p1) 
    return HugoniotVelocity(p2, p1, g);
  else
    return IsentropeVelocity(p2, p1, g);
}

double dHydrodynamicVelocity(double p2, double p1, double g)
{
  if (p2>p1)
    return dHugoniotVelocity(p2,p1,g);
  else
    return dIsentropeVelocity(p2,p1,g);
}

double PsAcoustic(double pl, double pr, double bl, double br, double g)
{
  return pl*pr*(2+(bl-br)*g/pow(g-1,1.5))/(pl+pr);
}

namespace {
/*! \brief Strong rarefaction approximation to the pressure at the interface of the Riemann problem
  \param pl Left pressure
  \param pr Right pressure
  \param bl Left dimensionless velocity
  \param br Right dimensionless velocity
  \param g Adiabatic index
  \return Pressure
 */
double PsStrongRarefaction(double pl, double pr,
			   double bl, double br,
			   double g)
{
  return pow(pow(pr,-2*sqrt(g-1)/g)*(1.+br)/(1.-br)+
	     pow(pl,-2*sqrt(g-1)/g)*(1.-bl)/(1.+bl),
	     -g/2/sqrt(g-1));
}
}

double LeftVelocity(double b0, double p, double p0, double g)
{
  return RelVelAdd(b0,-HydrodynamicVelocity(p,p0,g));
}

double dLeftVelocity(double b0, double p, double p0, double g)
{
  return dRelVelAdd(b0,-HydrodynamicVelocity(p,p0,g),
		    0, -dHydrodynamicVelocity(p,p0,g));
}

double RightVelocity(double b0, double p, double p0, double g)
{
  return RelVelAdd(b0,HydrodynamicVelocity(p,p0,g));
}

double dRightVelocity(double b0, double p, double p0, double g)
{
  return dRelVelAdd(b0,HydrodynamicVelocity(p,p0,g),
		    0,dHydrodynamicVelocity(p,p0,g));
}

double TranscendentalEquation(double bl, double br, 
			      double pl, double pr, 
			      double p, double g)
{
  return RightVelocity(br,p,pr,g) - LeftVelocity(bl,p,pl,g);
}

double dTranscendentalEquation(double bl, double br, double pl, double pr, double p, double g)
{
  return dRightVelocity(br,p,pr,g) - dLeftVelocity(bl,p,pl,g);
}

/*
  Class constructor
  bli, bri - Left and right dimensionless velocities
  pli, pri - Left and right pressures
  gi - Adiabatic index
 */
iTranscendentalEquation::iTranscendentalEquation(double bli, double bri,
						 double pli, double pri,
						 double gi):
  bl(bli),
  br(bri),
  pl(pli),
  pr(pri),
  g(gi) {}

double iTranscendentalEquation::operator()(double x) const
{
  return TranscendentalEquation(bl, br, pl, pr, x, g);
}

double iTranscendentalEquation::diff(double x) const
{
  return dTranscendentalEquation(bl, br, pl, pr, x, g);
}

double SolveTranscendentalEquation(double bl, double br, 
				   double pl, double pr, 
				   double g)
{
  iTranscendentalEquation ite(bl, br, pl, pr, g);
  double pg = PsAcoustic(pl, pr, bl, br, g);
  if (pg<0)
    pg = PsStrongRarefaction(pl, pr, bl, br, g);
  SVRelStep sc(1e-12);
  return NewtonRaphson(ite, pg, sc);
}

RiemannSolution RiemannSolve(double bl, double br,
			     double pl, double pr,
			     double g)
{
  double ps = SolveTranscendentalEquation(bl, br, pl, pr, g);
  double bs = 0.5*(LeftVelocity(bl,ps,pl,g)+
		   RightVelocity(br,ps,pr,g));
  RiemannSolution res;
  res.Pressure = ps;
  res.Velocity = bs;
  return res;
}
