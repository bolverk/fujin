#include <cassert>
#include <functional>
#include "hydrodynamic_variables.hpp"
#include "utilities.hpp"
#include "universal_error.hpp"

using std::function;

Conserved::Conserved(const array<double, 3>& source):
  array<double, 3>{source},
  Mass((*this)[0]),
  Momentum((*this)[1]),
  Energy((*this)[2]) {}

Conserved::Conserved(void):
  Conserved(array<double, 3>()) {}

Conserved::Conserved(double mass,
		     double momentum,
		     double energy):
  Conserved(array<double, 3>{mass, momentum, energy}) {}

Conserved& Conserved::operator=(const Conserved& source)
{
  array<double, 3>::operator=(source);
  return *this;
}

NewConserved::NewConserved(const array<double, 3>& source):
  array<double, 3>{source},
  mass{(*this)[0]},
  positive{(*this)[1]},
  negative{(*this)[2]} {}

NewConserved::NewConserved(void):
  NewConserved{array<double, 3>()} {}

NewConserved::NewConserved(double mass_i,
			   double positive_i,
			   double negative_i):
  NewConserved{array<double, 3>{mass_i, positive_i, negative_i}} {}

NewConserved::NewConserved(const NewConserved& source):
  NewConserved(array<double, 3>(source)){}

NewConserved& NewConserved::operator=(const NewConserved& source)
{
  array<double, 3>::operator=(source);
  return *this;
}

Primitive::Primitive(const array<double, 3>& source):
  array<double,3>(source),
  Density((*this)[0]),
  Pressure((*this)[1]),
  Celerity((*this)[2]) {}

Primitive::Primitive(void):
  Primitive(array<double, 3>()) {}

Primitive::Primitive(double density,
		     double pressure,
		     double celerity):
  Primitive(array<double, 3>{density,pressure,celerity}) {}

Primitive::Primitive(const Primitive& source):
  Primitive(array<double, 3>(source)) {}

Primitive& Primitive::operator=(const Primitive& source)
{
  array<double, 3>::operator=(source);
  return *this;
}

namespace {
  template<class T> T bin_op
  (const T& p1,
   const T& p2,
   function<double(double, double)> func)
  {
    T res;
    transform(p1.begin(),
	      p1.end(),
	      p2.begin(),
	      res.begin(),
	      func);
    return res;
  }

  template<class T> T unary_op
  (const T& p1,
   function<double(double)> func)
  {
    T res;
    transform(p1.begin(),
	      p1.end(),
	      res.begin(),
	      func);
    return res;
  }
}

Primitive operator+(Primitive const& p1,
		    Primitive const& p2)
{
  return bin_op(p1, p2, std::plus<double>());
}

Primitive operator-(Primitive const& p1,
		    Primitive const& p2)
{
  return bin_op(p1, p2, std::minus<double>());
}

Primitive operator/(Primitive const& p,
		    double d)
{
  return unary_op(p, [&d](double s){return s/d;});
}

Primitive operator*(double d,
		    Primitive const& p)
{
  return unary_op(p, [&d](double s){return s*d;});
}

HydroSnapshot::HydroSnapshot
(vector<double> const& edges_i,
 vector<Primitive> const& cells_i):
  edges(edges_i), cells(cells_i) {}

HydroSnapshot::HydroSnapshot(HydroSnapshot const& source):
  edges(source.edges),cells(source.cells) {}

HydroSnapshot& HydroSnapshot::operator=
(const HydroSnapshot& source)
{
  edges = source.edges;
  cells = source.cells;
  return *this;
}

HydroSnapshot operator+(HydroSnapshot const& hs1,
			HydroSnapshot const& hs2)
{
  return HydroSnapshot(hs1.edges+hs2.edges,
		       hs1.cells+hs2.cells);
}

HydroSnapshot operator-(HydroSnapshot const& hs1,
			HydroSnapshot const& hs2)
{
  return HydroSnapshot(hs1.edges-hs2.edges,
		       hs1.cells-hs2.cells);
}

HydroSnapshot operator*(double d,
			HydroSnapshot const& hs)
{
  return HydroSnapshot(d*hs.edges,d*hs.cells);
}

Conserved operator+(Conserved const& v1, 
		    Conserved const& v2)
{
  return bin_op(v1, v2, std::plus<double>());
}

NewConserved operator+(const NewConserved& v1, 
		       const NewConserved& v2)
{
  return bin_op(v1, v2, std::plus<double>());
}

Conserved operator-(Conserved const& v1, Conserved const& v2)
{
  return bin_op(v1, v2, std::minus<double>());
}

NewConserved operator-(const NewConserved& v1, const NewConserved& v2)
{
  return bin_op(v1, v2, std::minus<double>());
}

Conserved operator*(double s, Conserved const& v)
{
  return unary_op(v, [&s](double d){return d*s;});
}

NewConserved operator*(double s, const NewConserved& v)
{
  return unary_op(v, [&s](double d){return s*d;});
}

Conserved operator/(Conserved const& v, double s)
{
  return (1/s)*v;
}
