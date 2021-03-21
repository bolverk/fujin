#ifndef VAN_LEER_HPP
#define VAN_LEER_HPP 1

#include "spatial_reconstruction.hpp"

/*! \brief Spatial reconstruction based on Van Leer's method
  \details The description of the method is available in <a href="http://www.sciencedirect.com/science/article/pii/S0021999197957041">B. Van Leer, 'Towards the Ultimate Conservative Difference Scheme', J. Comp. Phys. 135, 229-248 (1977) </a>
 */
class VanLeer: public SpatialReconstruction
{
public:

  /*! \brief Class constructor
   */
  VanLeer(void);

  vector<std::pair<Primitive,Primitive> > interpolateAll
  (const NewHydroSnapshot<vector<double>, vector<Primitive> >& hs,
   double dt) const;
};

#endif // VAN_LEER_HPP
