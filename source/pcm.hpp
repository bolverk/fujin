#ifndef PCM_HPP
#define PCM_HPP 1

#include "utilities.hpp"

#include "spatial_reconstruction.hpp"

namespace pcm{
  template<template<class> class CE, template<class> class CP>
  class Interpolator: public Index2Member<std::pair<Primitive,Primitive> >
  {
  public:

    /*! \brief Class constructor
      \param hs HydroSnapshot
    */
    explicit Interpolator
    (const NewHydroSnapshot<CE, CP>& hs):
      hs_(hs) {}

    size_t getLength(void) const
    {
      return hs_.cells.size()+1;
    }

    std::pair<Primitive,Primitive> operator()(size_t i) const
    {
      if(i==0 || i==getLength()-1)
	return std::pair<Primitive,Primitive>();
      return std::pair<Primitive,Primitive>(hs_.cells[i-1],
					    hs_.cells[i]);
    }

  private:
    //! \brief Hydrodynamic snapshot
    const NewHydroSnapshot<CE, CP>& hs_;
  };
}

//! \brief Piecewise constant interpolation
template<template<class> class CE, template<class> class CP>
class PCM: public SpatialReconstruction<CE, CP>
{

public:

  //! \brief Class constructor
  PCM(void) {}

  CE<std::pair<Primitive,Primitive> > 
  interpolateAll
  (const NewHydroSnapshot<CE, CP>& hs,
   double /*dt*/) const override
  {
    return serial_generate<pair<Primitive, Primitive>, CE>
      (pcm::Interpolator<CE, CP>(hs));
  }
};

#endif // PCM_HPP
