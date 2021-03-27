#include <cassert>
#include "hydrodynamic_variables.hpp"
#include "utilities.hpp"
#include "universal_error.hpp"

using std::function;

Conserved::Conserved(const Conserved& source):
  Conserved(array<double, 3>(source)) {}

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

/*
HydroSnapshot::HydroSnapshot
(vector<double> const& edges_i,
 vector<Primitive> const& cells_i):
  pair<vector<double>, vector<Primitive> >{edges_i, cells_i},
  edges((*this).first),
  cells((*this).second) {}

HydroSnapshot::HydroSnapshot(HydroSnapshot const& source):
  HydroSnapshot(source.edges, source.cells) {}

HydroSnapshot& HydroSnapshot::operator=
(const HydroSnapshot& source)
{
  pair<vector<double>, vector<Primitive> >::operator=(source);
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
*/
