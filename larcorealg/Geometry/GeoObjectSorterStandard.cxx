////////////////////////////////////////////////////////////////////////
/// \file  GeoObjectSorterStandard.cxx
/// \brief Interface to algorithm class for sorting standard geo::XXXGeo objects
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////

#include "larcorealg/Geometry/GeoObjectSorterStandard.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/WireGeo.h"

#include <utility>

namespace {

  // Tolerance when comparing distances in geometry:
  constexpr double DistanceTol = 0.001; // cm

} // local namespace

namespace geo {

  //----------------------------------------------------------------------------
  // Define sort order for cryostats in standard configuration
  static bool sortAuxDetStandard(const AuxDetGeo& ad1, const AuxDetGeo& ad2)
  {
    // sort based off of GDML name, assuming ordering is encoded
    std::string const& ad1name = ad1.TotalVolume()->GetName();
    std::string const& ad2name = ad2.TotalVolume()->GetName();

    // assume volume name is "volAuxDet##"
    int ad1Num = atoi(ad1name.substr(9, ad1name.size()).c_str());
    int ad2Num = atoi(ad2name.substr(9, ad2name.size()).c_str());

    return ad1Num < ad2Num;
  }

  //----------------------------------------------------------------------------
  // Define sort order for cryostats in standard configuration
  static bool sortAuxDetSensitiveStandard(const AuxDetSensitiveGeo& ad1,
                                          const AuxDetSensitiveGeo& ad2)
  {
    // sort based off of GDML name, assuming ordering is encoded
    std::string ad1name = (ad1.TotalVolume())->GetName();
    std::string ad2name = (ad2.TotalVolume())->GetName();

    // assume volume name is "volAuxDetSensitive##"
    int ad1Num = atoi(ad1name.substr(9, ad1name.size()).c_str());
    int ad2Num = atoi(ad2name.substr(9, ad2name.size()).c_str());

    return ad1Num < ad2Num;
  }

  //----------------------------------------------------------------------------
  // Define sort order for cryostats in standard configuration
  static bool sortCryoStandard(const CryostatGeo& c1, const CryostatGeo& c2)
  {
    CryostatGeo::LocalPoint_t const local{0., 0., 0.};
    auto const xyz1 = c1.toWorldCoords(local);
    auto const xyz2 = c2.toWorldCoords(local);
    return xyz1.X() < xyz2.X();
  }

  //----------------------------------------------------------------------------
  // Define sort order for tpcs in standard configuration.
  static bool sortTPCStandard(const TPCGeo& t1, const TPCGeo& t2)
  {
    TPCGeo::LocalPoint_t const local{0., 0., 0.};
    auto const xyz1 = t1.toWorldCoords(local);
    auto const xyz2 = t2.toWorldCoords(local);

    // sort TPCs according to x
    return xyz1.X() < xyz2.X();
  }

  //----------------------------------------------------------------------------
  // Define sort order for planes in standard configuration
  static bool sortPlaneStandard(const PlaneGeo& p1, const PlaneGeo& p2)
  {
    PlaneGeo::LocalPoint_t const local{0., 0., 0.};
    auto const xyz1 = p1.toWorldCoords(local);
    auto const xyz2 = p2.toWorldCoords(local);

    // drift direction is negative, plane number increases in drift direction
    if (std::abs(xyz1.X() - xyz2.X()) > DistanceTol) return xyz1.X() > xyz2.X();

    //if same drift, sort by z
    if (std::abs(xyz1.Z() - xyz2.Z()) > DistanceTol) return xyz1.Z() < xyz2.Z();

    //if same z, sort by y
    return xyz1.Y() < xyz2.Y();
  }

  //----------------------------------------------------------------------------
  static bool sortWireStandard(WireGeo const& w1, WireGeo const& w2)
  {
    auto const [xyz1, xyz2] = std::make_pair(w1.GetCenter(), w2.GetCenter());

    //sort by z first
    if (std::abs(xyz1.Z() - xyz2.Z()) > DistanceTol) return xyz1.Z() < xyz2.Z();

    //if same z sort by y
    if (std::abs(xyz1.Y() - xyz2.Y()) > DistanceTol) return xyz1.Y() < xyz2.Y();

    //if same y sort by x
    return xyz1.X() < xyz2.X();
  }

  //----------------------------------------------------------------------------
  GeoObjectSorterStandard::GeoObjectSorterStandard(fhicl::ParameterSet const&) {}

  //----------------------------------------------------------------------------
  void GeoObjectSorterStandard::SortAuxDets(std::vector<geo::AuxDetGeo>& adgeo) const
  {
    std::sort(adgeo.begin(), adgeo.end(), sortAuxDetStandard);
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterStandard::SortAuxDetSensitive(
    std::vector<geo::AuxDetSensitiveGeo>& adsgeo) const
  {
    std::sort(adsgeo.begin(), adsgeo.end(), sortAuxDetSensitiveStandard);
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterStandard::SortCryostats(std::vector<geo::CryostatGeo>& cgeo) const
  {
    std::sort(cgeo.begin(), cgeo.end(), sortCryoStandard);
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterStandard::SortTPCs(std::vector<geo::TPCGeo>& tgeo) const
  {
    std::sort(tgeo.begin(), tgeo.end(), sortTPCStandard);
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterStandard::SortPlanes(std::vector<geo::PlaneGeo>& pgeo,
                                           geo::DriftDirection_t const driftDir) const
  {
    // sort the planes to increase in drift direction
    // The drift direction has to be set before this method is called.  It is set when
    // the CryostatGeo objects are sorted by the CryostatGeo::SortSubVolumes method
    if (driftDir == geo::kPosX)
      std::sort(pgeo.rbegin(), pgeo.rend(), sortPlaneStandard);
    else if (driftDir == geo::kNegX)
      std::sort(pgeo.begin(), pgeo.end(), sortPlaneStandard);
    else if (driftDir == geo::kUnknownDrift)
      throw cet::exception("TPCGeo") << "Drift direction is unknown, can't sort the planes\n";
  }

  //----------------------------------------------------------------------------
  void GeoObjectSorterStandard::SortWires(std::vector<geo::WireGeo>& wgeo) const
  {
    std::sort(wgeo.begin(), wgeo.end(), sortWireStandard);
  }

}
