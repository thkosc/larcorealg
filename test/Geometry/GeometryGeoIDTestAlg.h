/**
 * @file   GeometryGeoIDTestAlg.h
 * @brief  Tests the correct assignment of IDs to detector geometry objects
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   October 31, 2016
 */

#ifndef LARCORE_TEST_GEOMETRY_GEOMETRYGEOIDTESTALG_H
#define LARCORE_TEST_GEOMETRY_GEOMETRYGEOIDTESTALG_H

#include "larcorealg/Geometry/fwd.h"

namespace geo {

  class GeometryGeoIDTestAlg {
  public:
    explicit GeometryGeoIDTestAlg(geo::GeometryCore const* new_geo) : geom{new_geo} {}

    /// Executes the test
    unsigned int Run() const;

    /// @name All the ID iterator tests
    /// @{
    void CryostatGeoIDTest() const;
    void TPCGeoIDTest() const;
    void PlaneGeoIDTest() const;
    void WireGeoIDTest() const;
    /// @}

  private:
    GeometryCore const* geom = nullptr; ///< pointer to the geometry description

  }; // class GeometryGeoIDTestAlg

} // namespace geo

#endif // LARCORE_TEST_GEOMETRY_GEOMETRYGEOIDTESTALG_H
