////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorter.cxx
/// \brief Interface to algorithm class for sorting geo::XXXGeo objects
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/Geometry/OpDetGeo.h"

#include <algorithm>
#include <vector>

namespace {
  bool sortorderOpDets(const geo::OpDetGeo& t1, const geo::OpDetGeo& t2)
  {
    geo::OpDetGeo::LocalPoint_t const local{};
    auto const xyz1 = t1.toWorldCoords(local);
    auto const xyz2 = t2.toWorldCoords(local);

    if (xyz1.Z() != xyz2.Z()) return xyz1.Z() > xyz2.Z();
    if (xyz1.Y() != xyz2.Y()) return xyz1.Y() > xyz2.Y();
    return xyz1.X() > xyz2.X();
  }
}

namespace geo {
  void GeoObjectSorter::SortOpDets(std::vector<geo::OpDetGeo>& opdet) const
  {
    std::sort(opdet.begin(), opdet.end(), sortorderOpDets);
  }
}
