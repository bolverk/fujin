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

vector<double> distribute_vertices1
(const vector<double>& vertices)
{
#ifdef PARALLEL
  spdlog::debug("Inside parallel distribute vertices");
  const size_t cell_num = vertices.size()-1;
  const int rank = get_mpi_rank();
  const int size = get_mpi_size();
  spdlog::debug("cell_num {0}, rank {1}, size {2}",
		cell_num, rank, size);
    
  vector<int> partition(static_cast<size_t>(size), cell_num/size);
  for(size_t i=0;i<cell_num%size;++i)
    ++partition.at(i);
    
  vector<int> cumpar(partition.size()+1,0);
  for(size_t i=1;i<cumpar.size();++i)
    cumpar.at(i) = cumpar.at(i-1) + partition.at(i-1);

  const size_t low = cumpar.at(rank);
  const size_t high = cumpar.at(rank+1);

  vector<double> res(high-low+1,0);
  for(size_t i=0;i<res.size();++i)
    res.at(i) = vertices.at(low+i);
  spdlog::debug("low {0}, high {1}",
		res.front(),
		res.back());
  return res;
#else
  return vertices;
#endif // PARALLEL
}

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
