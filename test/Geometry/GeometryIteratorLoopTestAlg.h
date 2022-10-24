/**
 * @file   GeometryIteratorLoopTestAlg.h
 * @brief  Tests the correct iteration of the geo::Geometry iterators
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   May 7th, 2015
 */

#ifndef GEO_GEOMETRYITERATORLOOPTESTALG_H
#define GEO_GEOMETRYITERATORLOOPTESTALG_H

#include "larcorealg/Geometry/fwd.h"

namespace geo {

  class GeometryIteratorLoopTestAlg {
  public:
    GeometryIteratorLoopTestAlg(geo::GeometryCore const* new_geo) : geom{new_geo} {}

    unsigned int Run();

  private:
    GeometryCore const* geom = nullptr;
  };

} // namespace geo

#endif // GEO_GEOMETRYITERATORLOOPTESTALG_H
