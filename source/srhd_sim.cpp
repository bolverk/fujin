/*! \file srhd_sim.cpp
  \brief Relativistic, Godunov type, lagrangian 1d hydrodynamic simulation
  \author Almog Yalinewich
 */

#include <cmath>
#include <complex>
#include <cassert>
#include "srhd_sim.hpp"
#include "utilities.hpp"
#include "srhydro.hpp"
#include "universal_error.hpp"
#include "advanced_hydrodynamic_variables.hpp"
#ifdef PARALLEL
#include "mpi.h"
#endif // PARALLEL
#include "spdlog/spdlog.h"

using namespace std;

namespace{

#ifdef PARALLEL
  int get_mpi_rank(void)
  {
    int res;
    MPI_Comm_rank(MPI_COMM_WORLD, &res);
    return res;
  }

  int get_mpi_size(void)
  {
    int res;
    MPI_Comm_size(MPI_COMM_WORLD, & res);
    return res;
  }
#endif // PARALLEL
}
