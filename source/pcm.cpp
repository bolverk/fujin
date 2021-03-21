#include "pcm.hpp"
#include "advanced_hydrodynamic_variables.hpp"
#include "utilities.hpp"

PCM::PCM(void) {}

namespace {

  //! \brief Calculates the interpolated values
  class Interpolator: public Index2Member<std::pair<Primitive,Primitive> >
  {
  public:

    /*! \brief Class constructor
      \param hs HydroSnapshot
     */
    explicit Interpolator
    (const NewHydroSnapshot<vector<double>, vector<Primitive> >& hs):
      hs_(hs) {}

    size_t getLength(void) const
    {
      return hs_.cells.size()-1;
    }

    std::pair<Primitive,Primitive> operator()(size_t i) const
    {
      return std::pair<Primitive,Primitive>(hs_.cells[i],
					    hs_.cells[i+1]);
    }

  private:
    //! \brief Hydrodynamic snapshot
    const NewHydroSnapshot<vector<double>, vector<Primitive> >& hs_;
  };
}

vector<std::pair<Primitive,Primitive> > PCM::interpolateAll
(const NewHydroSnapshot<vector<double>, vector<Primitive> >& hs,
 double /*dt*/) const
{
  return serial_generate(Interpolator(hs));
}
