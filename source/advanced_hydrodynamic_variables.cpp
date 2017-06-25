#include "advanced_hydrodynamic_variables.hpp"
#include "utilities.hpp"

Conserved Primitive2Conserved(Primitive const& p, EquationOfState const& eos)
{
  // Lorentz factor
  const double lf = celerity2lorentz_factor(p.Celerity);
  // Enthalpy
  const double energy = eos.dp2e(p.Density, p.Pressure);
  const double h = (p.Pressure + energy)/p.Density;
  return Conserved(p.Density*lf,
		   h*p.Celerity,
		   h * lf - (p.Pressure/p.Density)/lf);
}

NewConserved primitive_to_new_conserved(const Primitive& p,
					const EquationOfState& eos)
{
  const double lf = celerity2lorentz_factor(p.Celerity);
  const double energy = eos.dp2e(p.Density, p.Pressure);
  const double h = (p.Pressure+energy)/p.Density;
  return NewConserved
    (p.Density*lf,
     p.Celerity > 0 ?
     h*(lf+p.Celerity)-(p.Pressure/p.Density)/lf :
     h/(lf-p.Celerity)-(p.Pressure/p.Density)/lf,
     p.Celerity < 0 ?
     h*(lf-p.Celerity)-(p.Pressure/p.Density)/lf :
     h/(lf+p.Celerity)-(p.Pressure/p.Density)/lf);
}

NewConserved old_to_new_conserved(const Conserved& c)
{
  return NewConserved(c.Mass,
		      c.Energy+c.Momentum,
		      c.Energy-c.Momentum);
}

Conserved new_to_old_conserved(const NewConserved& c)
{
  return Conserved(c.mass,
		   0.5*(c.positive - c.negative),
		   0.5*(c.positive + c.negative));
}

Conserved Primitive2Flux(Primitive const& hs, EquationOfState const& eos)
{
  const double energy = eos.dp2e(hs.Density, hs.Pressure);
  const double h = energy+hs.Pressure;
  const double lf = celerity2lorentz_factor(hs.Celerity);
  return Conserved(hs.Density*hs.Celerity,
		   (pow(hs.Celerity,2)*h+hs.Pressure)/(hs.Density*lf),
		   h*hs.Celerity/hs.Density);
}

NewConserved primitive_to_flux(const Primitive& hs,
			       const EquationOfState& eos)
{
  const double energy = eos.dp2e(hs.Density, hs.Pressure);
  const double h = energy + hs.Pressure;
  const double lf = celerity2lorentz_factor(hs.Celerity);
  const double w = hs.Celerity;
  const double d = hs.Density;
  const double p = hs.Pressure;
  return NewConserved(hs.Density*hs.Celerity,
		      w > 0 ?
		      (w*h/d)*(lf+w)/lf+(p/d)/lf :
		      (w*h/d)/(lf-w)/lf+(p/d)/lf,
		      w < 0 ?
		      (w*h/d)*(lf-w)/lf-(p/d)/lf :
		      (w*h/d)/(lf+w)/lf-(p/d)/lf);
}

Conserved Primitive2Conserved_pv(Primitive const& p,
				 const EquationOfState& eos)
{
  Conserved res = Primitive2Conserved(p, eos);
  res.Momentum *= res.Mass;
  res.Energy *= res.Mass;
  return res;
}

NewConserved primitive_to_conserved_pv(const Primitive& p,
				       const EquationOfState& eos)
{
  NewConserved res = primitive_to_new_conserved(p,eos);
  res.positive *= res.mass;
  res.negative *= res.mass;
  return res;
}

Conserved Primitive2Flux_pv(Primitive const& p,
			    const EquationOfState& eos)
{
  const double mass = p.Density*celerity2lorentz_factor(p.Celerity);
  Conserved res = Primitive2Flux(p,eos);
  res.Momentum *= mass;
  res.Energy *= mass;
  return res;
}

NewConserved primitive_to_flux_pv(const Primitive& hs,
				  const EquationOfState& eos)
{
  const double mass = hs.Density*celerity2lorentz_factor(hs.Celerity);
  NewConserved res = primitive_to_flux(hs,eos);
  res.positive *= mass;
  res.negative *= mass;
  return res;
}
