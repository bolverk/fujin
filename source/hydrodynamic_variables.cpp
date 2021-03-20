#include <cassert>
#include "hydrodynamic_variables.hpp"
#include "utilities.hpp"
#include "universal_error.hpp"

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

Primitive operator+(Primitive const& p1,
		    Primitive const& p2)
{
  return Primitive(p1.Density+p2.Density,
		   p1.Pressure+p2.Pressure,
		   p1.Celerity+p2.Celerity);
}

Primitive operator-(Primitive const& p1,
		    Primitive const& p2)
{
  return Primitive(p1.Density-p2.Density,
		   p1.Pressure-p2.Pressure,
		   p1.Celerity-p2.Celerity);
}

Primitive operator/(Primitive const& p,
		    double d)
{
  return Primitive(p.Density/d,
		   p.Pressure/d,
		   p.Celerity/d);
}

Primitive operator*(double d,
		    Primitive const& p)
{
  return Primitive(d*p.Density,
		   d*p.Pressure,
		   d*p.Celerity);
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
  return Conserved(v1.Mass+v2.Mass,
		   v1.Momentum+v2.Momentum,
		   v1.Energy+v2.Energy);
}

NewConserved operator+(const NewConserved& v1, 
		       const NewConserved& v2)
{
  return NewConserved(v1.mass+v2.mass,
		      v1.positive+v2.positive,
		      v1.negative+v2.negative);
}

Conserved operator-(Conserved const& v1, Conserved const& v2)
{
  return Conserved(v1.Mass-v2.Mass,
		   v1.Momentum-v2.Momentum,
		   v1.Energy-v2.Energy);
}

NewConserved operator-(const NewConserved& v1, const NewConserved& v2)
{
  return NewConserved(v1.mass-v2.mass,
		      v1.positive-v2.positive,
		      v1.negative-v2.negative);
}

Conserved operator*(double s, Conserved const& v)
{
  return Conserved(s*v.Mass,
		   s*v.Momentum,
		   s*v.Energy);
}

NewConserved operator*(double s, const NewConserved& v)
{
  return NewConserved(s*v.mass,
		      s*v.positive,
		      s*v.negative);
}

Conserved operator/(Conserved const& v, double s)
{
  return (1/s)*v;
}
