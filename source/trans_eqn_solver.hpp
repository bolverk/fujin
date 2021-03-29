/*! \file trans_eqn_solver.hpp
  \brief Various solvers for single variable transcendental equations
  \author Almog Yalinewich
 */
#ifndef TRANS_EQN_SOLVER
#define TRANS_EQN_SOLVER 1

/*! \brief A single variable function and its derivative
  \details Input for the Newton raphson method
 */
class SVFunction
{
public:
  /*! \brief Evaluates the function
    \param x Independent variable    
    \return Double
   */
  virtual double operator()(double x) const = 0;

  virtual ~SVFunction(void);
};

//! \brief Single variable differentiable function
class SVDifferentiable: public SVFunction
{
public:

  virtual double operator()(double x) const override = 0;

  /*! \brief Evaluates the derivative
    \param x Argument
    \return derivative at x
   */
  virtual double diff(double x) const = 0;

  virtual ~SVDifferentiable(void) override;
};

/*! \brief Stopping condition for the solvers
 */
class SVStopCond
{
public:

  /*! \brief Determines whether the solution meets a certain presicion criterion
    \param x Current value of the independent varaible
    \param f Current value of the function
    \param svf Pointer to the function
    \return True if the convergence criterion is met, false otherwise
   */
  virtual bool operator()(double x, double f, SVFunction const& svf) const = 0;

  virtual ~SVStopCond(void);
};

/*! \brief Hybrid bisection / Newton - Raphson method
  \param svf Transcendental equation
  \param xl Left (lower) bracket
  \param xr Right (upper) bracket
  \param sc Stopping condition
  \return Value of the function
 */
double NRSafe(SVDifferentiable const& svf, 
	      double xl, double xr, 
	      const SVStopCond& sc);

/*! \brief Bisection method
  \param svf Transcendental equation
  \param xl Left (lower) bracket
  \param xr Right (upper) bracket
  \param sc Stopping condition
  \return solution to the equation
*/
double Bisection(SVFunction const& svf, double xl, double xr, 
		 const SVStopCond& sc);

/*! \brief Newton Raphson method
  \param svf Transcendental equation
  \param xg First guess of the solution
  \param sc Stopping condition
  \return Solution to the equation
 */
double NewtonRaphson(SVDifferentiable const& svf, double xg, 
		     const SVStopCond& sc);

//! \brief Stopping condition according to ratio between step size magnitude of the solution
class SVRelStep: public SVStopCond
{
public:

  /*! \brief Class constructor
    \param tol Tolerance
   */
  explicit SVRelStep(double tol);

  bool operator()(double x, double f, SVFunction const& svf) const override;

  /*! \brief Returns the number of iterations
   */
  int get_iter(void) const;

private:

  //! \brief Tolerance
  const double tol_;

  //! \brief Number of iterations
  mutable int iter_;

  //! \brief Value of independent variable in the previous cycle
  mutable double xold_;
};

#endif // TRANS_EQN_SOLVER
