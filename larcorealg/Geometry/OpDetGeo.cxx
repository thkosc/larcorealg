////////////////////////////////////////////////////////////////////////
/// \file  OpDetGeo.cxx
/// \brief Encapsulate the geometry of an OpDet
///
/// \author  bjpjones@mit.edu
////////////////////////////////////////////////////////////////////////

// class header
#include "larcorealg/Geometry/OpDetGeo.h"

// LArSoft libraries
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect::makeFromCoords()
#include "larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h" // util::pi()

// ROOT libraries
#include "TGeoManager.h"
#include "TGeoTube.h"
#include "TGeoNode.h"

// C/C++ standard libraries
#include <cmath>


namespace {
  template <typename T>
  inline T sqr(T v) { return v*v; }
} // local namespace

namespace geo{

  //-----------------------------------------
  OpDetGeo::OpDetGeo(std::vector<const TGeoNode*>& path, int depth)
    : fTrans(path, depth)
  {
    fOpDetNode = path[depth];
    
    fCenter = toWorldCoords(geo::origin<LocalPoint_t>());
    
  }

  //......................................................................

  /// Return the center position of an opdet
  /// \param xyz : 3-D array. The returned location.
  /// \param localz : Distance along the length of the volume
  /// (cm). Default is center of wire
  void OpDetGeo::GetCenter(double* xyz, double localz) const
  {
    double xyzLocal[3] = {0.,0.,localz};
    this->LocalToWorld(xyzLocal, xyz);
  }

  //......................................................................

  double OpDetGeo::RMax() const 
  {
    return ((TGeoTube*)fOpDetNode->GetVolume()->GetShape())->GetRmax();
  }
    
  //......................................................................

  double OpDetGeo::HalfL() const 
  {
    return ((TGeoTube*)fOpDetNode->GetVolume()->GetShape())->GetDZ();
  }

  //......................................................................

  double OpDetGeo::RMin() const 
  {
    return ((TGeoTube*)fOpDetNode->GetVolume()->GetShape())->GetRmin();
  }

  //......................................................................
  double OpDetGeo::ThetaZ() const
  {
    auto const& center = GetCenter();
    auto const& end = toWorldCoords(LocalPoint_t{ 0.0, 0.0, HalfL() });
    
    // TODO change this into something generic
    //either y or x will be 0, so ading both will always catch the right
    //one
    double angle = (end.Y()-center.Y()+end.X()-center.X()) / 
      std::abs(end.Y()-center.Y()+center.X()-end.X()) * 
      std::acos((end.Z() - center.Z())/HalfL());
    if (angle < 0) angle += util::pi(); 
    return angle;
  }

  //......................................................................
  double OpDetGeo::ThetaZ(bool degree) const
    { return degree? util::RadiansToDegrees(ThetaZ()): ThetaZ(); }
  
  

  //......................................................................
  double OpDetGeo::DistanceToPoint(geo::Point_t const& point) const
    { return (point - GetCenter()).R(); }
  double OpDetGeo::DistanceToPoint(double const* xyz) const
    { return DistanceToPoint(geo::vect::makeFromCoords<geo::Point_t>(xyz)); }


  //......................................................................
  double OpDetGeo::CosThetaFromNormal(geo::Point_t const& point) const {
    auto const& local = toLocalCoords(point);
    return local.Z() / local.R();
  }
  double OpDetGeo::CosThetaFromNormal(double const* xyz) const
    { return CosThetaFromNormal(geo::vect::makeFromCoords<geo::Point_t>(xyz)); }


}
////////////////////////////////////////////////////////////////////////
