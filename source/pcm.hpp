#ifndef PCM_HPP
#define PCM_HPP 1

#include "spatial_reconstruction.hpp"

template<class T> using simple_vector = vector<T>;

//! \brief Piecewise constant interpolation
class PCM: public SpatialReconstruction<simple_vector, simple_vector>
{

public:

  //! \brief Class constructor
  PCM(void);

  vector<std::pair<Primitive,Primitive> > 
  interpolateAll
  (const NewHydroSnapshot<simple_vector, simple_vector>& hs,
   double dt) const;
};

#endif // PCM_HPP
