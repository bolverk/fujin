/*! \file spatial_distribution.hpp
  \brief Abstract class for initial conditions
  \author Almog Yalinewich
 */

#include "spatial_distribution.hpp"
#include <cmath>

SpatialDistribution::~SpatialDistribution(void) {}

Uniform::Uniform(double iValue):
  Value(iValue) {}

double Uniform::operator()(double /*x*/) const
{
  return Value;
}

Step::Step(double iValue1, double iValue2, 
	   double iStepPosition):
  Value1(iValue1),
  Value2(iValue2),
  StepPosition(iStepPosition) {}

double Step::operator()(double x) const
{
  if (x>StepPosition)
    return Value2;
  else
    return Value1;
}

TwoSteps::TwoSteps(double ivl, double ilip,
		   double ivm, double irip,
		   double ivr):
  vl(ivl),
  lip(ilip),
  vm(ivm),
  rip(irip),
  vr(ivr) {}

double TwoSteps::operator()(double x) const
{
  if (x<lip)
    return vl;
  else if (x>rip)
    return vr;
  else
    return vm;
}

SineWave::SineWave(double amplitude,
		   double wavelength,
		   double phase,
		   double offset):
  amp_(amplitude),
  k_(2*M_PI/wavelength),
  ph_(phase),
  offset_(offset) {}

double SineWave::operator()(double x) const
{
  return amp_*sin(k_*x+ph_)+offset_;
}
