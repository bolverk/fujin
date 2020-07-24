/*
  Verifies calculation of enthalpy
 */

#include <fenv.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <boost/foreach.hpp>
#include "ideal_gas.hpp"
#include "hydrodynamic_variables.hpp"
#include "advanced_hydrodynamic_variables.hpp"
#include "utilities.hpp"
#include "hdf5_utils.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

using namespace std;

namespace {
  vector<double> w_range(void)
  {
    vector<double> res;
    for(int i=-6;i<7;++i)
      res.push_back(-pow(10,-static_cast<double>(i)));
    res.push_back(0);
    for(int i=-6;i<7;++i)
      res.push_back(pow(10,static_cast<double>(i)));
    return res;
  }

  vector<double> p_range(void)
  {
    vector<double> res;
    for(int i=-6;i<7;++i)
      res.push_back(pow(10,static_cast<double>(i)));
    return res;
  }
}

int main()
{
#ifdef PARALLEL
  MPI_Init(NULL, NULL);
#endif // PARALLEL
  feenableexcept(FE_INVALID   | 
		 FE_DIVBYZERO | 
		 FE_OVERFLOW  | 
		 FE_UNDERFLOW);

  const double g = 4./3.;
  const IdealGas eos(4./3.);

  vector<double> pressure;
  vector<double> celerity;
  vector<double> reconstructed;
  BOOST_FOREACH(double p, p_range())
    {
      BOOST_FOREACH(double w, w_range())
	{
	  const Primitive hs(1,p,w);
	  const NewConserved cv = primitive_to_new_conserved(hs,eos);
	  const double dh = conserved2enthalpy_diff
	    (cv.positive, cv.negative, g);
	  pressure.push_back(p);
	  celerity.push_back(w);
	  reconstructed.push_back(dh*(1.-1./g));
	}
    }

#ifdef PARALLEL
  if(get_mpi_rank()==0)
#endif // PARALLEL
  HDF5Shortcut("results.h5")
    ("pressure",pressure)
    ("celerity",celerity)
    ("reconstructed",reconstructed);

  // Finalise
  ofstream("test_terminated_normally.res").close();

#ifdef PARALLEL
  MPI_Finalize();
#endif // PARALLEL
  return 0;
}
