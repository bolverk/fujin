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

  double calcArea(double radius) const override;

  double calcVolume(double radius) const override;
};

//! \brief Cylindrical geometry
class Cylindrical: public Geometry
{
public:

  Cylindrical(void);

  double calcArea(double radius) const override;

  double calcVolume(double radius) const override;
};

//! \brief Spherical geometry
class Spherical: public Geometry
{
public:

  Spherical(void);

  double calcArea(double radius) const override;

  double calcVolume(double radius) const override;
};

#endif // GEOMETRY_HPP
