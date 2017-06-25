#ifndef PCM_HPP
#define PCM_HPP 1

#include "spatial_reconstruction.hpp"

//! \brief Piecewise constant interpolation
class PCM: public SpatialReconstruction
{

public:

  //! \brief Class constructor
  PCM(void);

  vector<std::pair<Primitive,Primitive> > 
  interpolateAll
  (const HydroSnapshot& hs,
   double dt) const;
};

#endif // PCM_HPP
