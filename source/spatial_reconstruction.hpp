/*! \file spatial_reconstruction.hpp
  \brief Interpolation methods for the hydrodynamic states near the interface
  \author Almog Yalinewich
 */

#ifndef SPATIAL_RECONSTRUCTION_HPP
#define SPATIAL_RECONSTRUCTION_HPP

#include <vector>
#include "hydrodynamic_variables.hpp"

using std::min;

//! \brief Base class for interpolations
class SpatialReconstruction
{
public:

  /*! \brief Interpolate all cell values
    \param hs Hydrodynamic snapshot
    \param dt Time step
    \return List of interpolated values on both sides of each interface
   */
  virtual vector<std::pair<Primitive,Primitive> > 
  interpolateAll
  (const HydroSnapshot& hs,
   double dt) const = 0;

  virtual ~SpatialReconstruction(void);
};

#endif
