////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorter.h
/// \brief Interface to algorithm class for sorting geo::XXXGeo objects
/// \ingroup Geometry
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_AUXDETGEOOBJECTSORTER_H
#define GEO_AUXDETGEOOBJECTSORTER_H

#include <vector>

#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

#include "fhiclcpp/ParameterSet.h"

namespace geo{

  class AuxDetGeo;
  class AuxDetSensitiveGeo;

  /// \ingroup Geometry
  class AuxDetGeoObjectSorter {

  public:

    virtual ~AuxDetGeoObjectSorter() = default;

    virtual void SortAuxDets        (std::vector<geo::AuxDetGeo>          & adgeo)   const = 0;
    virtual void SortAuxDetSensitive(std::vector<geo::AuxDetSensitiveGeo> & adsgeo)  const = 0;

  };

}

#endif // GEO_GEOOBJECTSORTER_H
