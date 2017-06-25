#include <cmath>
#include <cassert>
#include "trans_eqn_solver.hpp"
#include "utilities.hpp"

using namespace std;

//! \brief An interval that contains a root of an equation
class Brackets
{
public:

  /*! \brief class constructor
    \param xl Left bracket
    \param fl Value of the function at xl
    \param xr Right bracket
    \param fr Value of the function at xr
  */
  Brackets(double xl, double fl, double xr, double fr):
    xl_(xl), fl_(fl), xr_(xr), fr_(fr)
  {
    assert(fl*fr<0);
  }

  /*! \brief Divides the space inside the bracket, and retains partition that contains the root
    \param xm x value inside he brackets
    \param fm Value of the function at xm
  */
  void Divide(double xm, double fm)
  {
    if(fm*fr_>0){
      xr_ = xm;
      fr_ = fm;
    }
    else if(fm*fl_>0){
      xl_ = xm;
      fl_ = fm;
    }
    else{
      assert(2<1 && "Unresolved boundary");
    }
  }

  /*! \brief Checks whether a value is inside the brackets
    \param xm x Value
    \return True if within bounds, false otherwise
  */
  bool IsBracketed(double xm) const
  {
    return (xm>=xl_&&xm<=xr_);
  }

  //! Width of the bracket
  double width(void) const
  {
    return xr_ - xl_;
  }

  //! \brief Position of left bracket
  double get_xl(void) const
  {
    return xl_;
  }

  //! \brief Position of right bracket
  double get_xr(void) const
  {
    return xr_;
  }

  //! \brief Value of the function at the left bracket
  double get_fl(void) const
  {
    return fl_;
  }

  //! \brief Value of the function at the right bracket
  double get_fr(void) const
  {
    return fr_;
  }

private:

  //! \brief Position of the left bracket
  double xl_;

  //! \brief Value of the function at the left bracket
  double fl_;

  //! \brief Position of the right bracket
  double xr_;

  //! \brief Value of the function at the right bracket
  double fr_;
};

namespace {

  /*! \brief Returns the middle of bracketed zone
    \param b Bracketed zone
    \return Position of the middle
   */
  double middle(Brackets const& b)
  {
    const double xr = b.get_xr();
    const double xl = b.get_xl();
    return 0.5*(xr+xl);
  }
}

SVFunction::~SVFunction(void) {}

SVDifferentiable::~SVDifferentiable(void) {}

SVStopCond::~SVStopCond(void) {}

double NRSafe(SVDifferentiable const& svf, 
	      double xl, double xr, 
	      const SVStopCond& sc)
{
  // Initialisation
  const double atol = 1e-15;
  const int max_iter = 1000;
  int iter = 0;

  // Main process
  double fl = svf(xl);
  double fr = svf(xr);
  if(abs(fl)<abs(atol*fr))
    return xl;
  if(abs(fr)<abs(atol*fl))
    return xr;

  Brackets brackets(xl, fl, xr, fr);
  double xm = middle(brackets);
  double fm = svf(xm);
  double dfm = svf.diff(xm);
  double dx_old = brackets.width();
  double dx = dx_old;
  while(!sc(xm,fm,svf)){
    if(!brackets.IsBracketed(xm-fm/dfm)|| // If not bracketed
       abs(2*fm)>abs(dx_old*dfm)){ // Or slow convergence
      dx_old = dx;
      dx = brackets.width()/2;
      xm = middle(brackets);
    }
    else{
      dx_old = dx;
      dx = -fm/dfm;
      xm += dx;
    }
    fm = svf(xm);
    dfm = svf.diff(xm);
    if(effectively_zero(fm))
      return xm;
    else
      brackets.Divide(xm,fm);

    iter++;
    assert(iter<max_iter);
  }
  return xm;
}

double Bisection(SVFunction const& svf, double xl, double xr, 
		 const SVStopCond& sc)
{
  // Initialisation
  const int max_iter = 100;
  int iter = 0;

  // Main process
  Brackets brackets(xl, svf(xl),
		    xr, svf(xr));
  double xm = middle(brackets);
  double fm = svf(xm);
  while(!sc(xm,fm,svf)){
    brackets.Divide(xm,fm);
    xm = middle(brackets);
    fm = svf(xm);
    iter++;
    assert(iter<max_iter);
  }
  return xm;
}

double NewtonRaphson(const SVDifferentiable& svf, double xg, 
		     const SVStopCond& sc)
{
  // Initialisation
  const int max_iter = 100;
  int iter = 0;
  const double xtol = 1e-15;

  // Main process
  double x = xg;
  double f = svf(x);
  
  while(!sc(x,f,svf)){
    f = svf(x);
    double df = svf.diff(x);
    if(abs(f)<xtol)
      return x;
    double dx = f/df;
    x -= dx;

    iter++;
    assert(iter<max_iter);
  }
  return x;
}

SVRelStep::SVRelStep(double tol):
  tol_(tol),iter_(0),xold_(0) {}

bool SVRelStep::operator()
(double x, double /*f*/, SVFunction const& /*svf*/) const
{
  if(2>iter_){
    iter_++;
    xold_ = x;
    return false;
  }
  else{
    double dx = x-xold_;
    iter_++;
    xold_ = x;
    return abs(dx/x)<tol_;
  }    
}

int SVRelStep::get_iter(void) const
{
  return iter_;
}
