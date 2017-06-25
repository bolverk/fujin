/*! \file spatial_distribution.hpp
  \brief Abstract class for initial conditions
  \author Almog Yalinewich
 */

#ifndef SPATIAL_DISTRIBUTION_HPP
#define SPATIAL_DISTRIBUTION_HPP 1

//! \brief Base class for initial conditions
class SpatialDistribution
{
public:
  /*! \brief Calculates initial conditions
    \param x Position
    \result Value of the function at x
   */
  virtual double operator()(double x) const = 0;

  virtual ~SpatialDistribution(void);
};

//! \brief Uniform distribution
class Uniform: public SpatialDistribution
{
private:
  //! \brief Value at each point
  double Value;
public:
  /*! \brief Class constructor
    \param iValue Value at each point
   */
  Uniform(double iValue);

  /*! \brief Calculates the value at each point
    \param x Position
    \return Value at x
   */
  double operator()(double x) const;
};

//! \brief Step distribution
class Step: public SpatialDistribution
{
private:
  //! \brief Value befor the step
  double Value1;
  //! \brief Value after the step
  double Value2;
  //! \brief Step position
  double StepPosition;
public:
  /*! \brief Class constructor
    \param iValue1 Value befor the step
    \param iValue2 Value after the step
    \param iStepPosition Step position
   */
  Step(double iValue1, double iValue2, 
       double iStepPosition);
  
  /*! \brief Evaluates the distribution at a certain spot
    \param x Position
    \return Value of the distribution at x
   */
  double operator()(double x) const;
};

//! \brief Two steps
class TwoSteps: public SpatialDistribution
{
private:
  //! \brief Value of the left step
  double vl;
  //! \brief Positions of the left interface
  double lip;
  //! \brief Value of the middle step
  double vm;
  //! \brief Position of the right interface
  double rip;
  //! \brief Value of the right step
  double vr;
public:
  /*! \brief Class constructor
    \param ivl Value on the left step
    \param ilip Position of the left interface
    \param ivm Value on the middle step
    \param irip Position of the right interface
    \param ivr Value on the right step
   */
  TwoSteps(double ivl, double ilip,
	   double ivm, double irip,
	   double ivr);
  /*! \brief Evaluates the distribution at a certain spot
    \param x Position
    \return Value of the distribution at x
   */
  double operator()(double x) const;  
};

//! \brief Sine wave distribution
class SineWave: public SpatialDistribution
{
public:

  /*! \brief Class constructor
    \param amplitude Amplitude
    \param wavelength Wavelength
    \param phase Phase
    \param offset Offset
   */
  SineWave(double amplitude,
	   double wavelength,
	   double phase,
	   double offset);

  double operator()(double x) const;

private:

  //! \brief Amplitude
  double amp_;

  //! \brief Wavenumber
  double k_;

  //! \brief Phase
  double ph_;

  //! \brief Offset
  double offset_;
};

#endif
