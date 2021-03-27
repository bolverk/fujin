#include <cassert>
#include "periodic.hpp"
#ifdef PARALLEL
#include "parallel_helper.hpp"
#endif // PARALLEL

//Periodic::Periodic(RiemannSolver const& rs):
//  rs_(rs) {}

#ifdef PARALLEL

namespace {

  vector<double> serialise_primitive
  (const Primitive& p)
  {
    vector<double> res(3,0);
    res.at(0) = p.Density;
    res.at(1) = p.Pressure;
    res.at(2) = p.Celerity;
    return res;
  }
  
  Primitive unserialise_primitive
    (const vector<double>& v)
  {
    assert(v.size()==3);
    Primitive res;
    res.Density = v.at(0);
    res.Pressure = v.at(1);
    res.Celerity = v.at(2);
    return res;
  }
}

RiemannSolution Periodic::operator()
  (bool idx,
   vector<Primitive> const& cells) const
{
  if(get_mpi_rank()==0){
    MPI_Request req;
    MPI_Isend(&serialise_primitive(cells.front()).front(),
	      3,
	      MPI_DOUBLE,
	      get_mpi_size()-1,
	      0,
	      MPI_COMM_WORLD,
	      &req);
    vector<double> buffer(3,0);
    MPI_Recv(&buffer.front(),
	     3,
	     MPI_DOUBLE,
	     get_mpi_size()-1,
	     0,
	     MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
    return rs_(unserialise_primitive(buffer),
	       cells.front());
  }
  if(get_mpi_rank()==get_mpi_size()-1){
    MPI_Request req;
    MPI_Isend(&serialise_primitive(cells.back()).front(),
	      3,
	      MPI_DOUBLE,
	      0,
	      0,
	      MPI_COMM_WORLD,
	      &req);
    vector<double> buffer(3,0);
    MPI_Recv(&buffer.front(),
	     3,
	     MPI_DOUBLE,
	     0,
	     0,
	     MPI_COMM_WORLD,
	     MPI_STATUS_IGNORE);
    MPI_Wait(&req, MPI_STATUS_IGNORE);
    return rs_(cells.back(),
	       unserialise_primitive(buffer));
  }
  throw("something bad happened in period boundary conditions");
}

#else

/*
RiemannSolution Periodic::operator()
(bool side,
 vector<Primitive> const& cells) const
{
  return side ? rs_(cells.front(),cells.back()) :
    rs_(cells.back(),cells.front());
    
}
*/

#endif // PARALLEL
