#include <cmath>
#include "geometry.hpp"

Geometry::~Geometry(void) {}

Planar::Planar(void) {}

double Planar::calcArea(double /*radius*/) const
{
  return 1;
}

double Planar::calcVolume(double radius) const
{
  return radius;
}

Cylindrical::Cylindrical(void) {}

double Cylindrical::calcArea(double radius) const
{
  return 2*M_PI*radius;
}

double Cylindrical::calcVolume(double radius) const
{
  return M_PI*pow(radius,2);
}

Spherical::Spherical(void) {}

double Spherical::calcArea(double radius) const
{
  return 4*M_PI*pow(radius,2);
}

double Spherical::calcVolume(double radius) const
{
  return (4.*M_PI/3.)*pow(radius,3);
}
