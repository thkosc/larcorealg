/**
 * @file   ChannelMapStandardTestAlg.h
 * @brief  Tests the standard channel mapping algorithm.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   June 26th, 2015
 */

#ifndef GEO_CHANNELMAPSTANDARDTESTALG_H
#define GEO_CHANNELMAPSTANDARDTESTALG_H

#include "larcorealg/Geometry/fwd.h"

namespace geo {

  class ChannelMapStandardTestAlg {
  public:
    explicit ChannelMapStandardTestAlg(geo::GeometryCore const* new_geo) : geom{new_geo} {}

    /// Executes the test
    unsigned int Run();

    void TPCsetMappingTest() const;
    void ROPMappingTest() const;
    void ChannelMappingTest() const;

  private:
    GeometryCore const* geom = nullptr;

  }; // class ChannelMapStandardTestAlg

} // namespace geo

#endif // GEO_CHANNELMAPSTANDARDTESTALG_H
