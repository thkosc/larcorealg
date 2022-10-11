////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterStandard.h
/// \brief Interface to algorithm class for standard sorting of geo::XXXGeo objects
/// \ingroup Geometry
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_GEOOBJECTSORTERSTANDARD_H
#define GEO_GEOOBJECTSORTERSTANDARD_H

#include <vector>

#include "fhiclcpp/fwd.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h" // for DriftDirec...

namespace geo {

  class AuxDetGeo;
  class AuxDetSensitiveGeo;
  class CryostatGeo;
  class PlaneGeo;
  class TPCGeo;
  class WireGeo;

  /// @ingroup Geometry
  class GeoObjectSorterStandard : public GeoObjectSorter {
  public:
    GeoObjectSorterStandard(fhicl::ParameterSet const& p);

    void SortAuxDets(std::vector<geo::AuxDetGeo>& adgeo) const override;
    void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo>& adsgeo) const override;
    void SortCryostats(std::vector<geo::CryostatGeo>& cgeo) const override;
    void SortTPCs(std::vector<geo::TPCGeo>& tgeo) const override;
    void SortPlanes(std::vector<geo::PlaneGeo>& pgeo,
                    geo::DriftDirection_t driftDir) const override;
    void SortWires(std::vector<geo::WireGeo>& wgeo) const override;
  };

}

#endif // GEO_GEOOBJECTSORTERSTANDARD_H
