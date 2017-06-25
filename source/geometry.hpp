#ifndef GEOMETRY_HPP
#define GEOMETRY_HPP 1

//! \brief Defines the geometry
class Geometry
{
public:

  /*! \brief Calculates the area
    \param radius Radius
    \return Area
   */
  virtual double calcArea(double radius) const = 0;

  /*! \brief Calculates the volume
    \param radius Radius
    \return Volume
   */
  virtual double calcVolume(double radius) const = 0;

  virtual ~Geometry(void);
};

//! \brief Planar geometry
class Planar: public Geometry
{
public:

  Planar(void);

  double calcArea(double radius) const;

  double calcVolume(double radius) const;
};

//! \brief Cylindrical geometry
class Cylindrical: public Geometry
{
public:

  Cylindrical(void);

  double calcArea(double radius) const;

  double calcVolume(double radius) const;
};

//! \brief Spherical geometry
class Spherical: public Geometry
{
public:

  Spherical(void);

  double calcArea(double radius) const;

  double calcVolume(double radius) const;
};

#endif // GEOMETRY_HPP
