#ifdef PARALLEL

#ifndef PARALLEL_HELPER
#define PARALLEL_HELPER 1

#include "mpi.h"

int get_mpi_rank(void);

int get_mpi_size(void);

#endif // PARALLEL_HELPER

#endif // PARALLEL
