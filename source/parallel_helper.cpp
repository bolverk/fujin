#ifdef PARALLEL

#include "parallel_helper.hpp"

int get_mpi_rank(void)
{
  int res = 0;
  MPI_Comm_rank(MPI_COMM_WORLD, &res);
  return res;
}

int get_mpi_size(void)
{
  int res = 0;
  MPI_Comm_size(MPI_COMM_WORLD, &res);
  return res;
}


#endif // PARALLEL
