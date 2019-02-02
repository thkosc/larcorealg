////////////////////////////////////////////////////////////////////////
/// \file  larcorealg/Geometry/OpDetGeo.h
/// \brief Encapsulate the geometry of an optical detector
/// \ingroup Geometry
///
/// \author  bjpjones@mit.gov
////////////////////////////////////////////////////////////////////////

#ifndef LARCOREALG_GEOMETRY_OPDETGEO_H
#define LARCOREALG_GEOMETRY_OPDETGEO_H

// LArSoft libraries
#include "larcorealg/Geometry/TransformationMatrix.h"
#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// ROOT libraries
#include "TGeoMatrix.h" // TGeoHMatrix
#include "TGeoTube.h"
#include "TGeoBBox.h"
#include "TClass.h"

// C/C++ standard libraries
#include <vector>
#include <string>
#include <array>


// forward declarations
class TGeoNode;

namespace geo {

  /// @ingroup Geometry
  class OpDetGeo {
  public:
    
    /// @{
    /**
     * @name Types for geometry-local reference vectors.
     * 
     * These types represents points and displacement vectors in the reference
     * frame defined in the optical detector geometry box from the GDML geometry
     * description.
     * 
     * No alias is explicitly defined for the LArSoft global vector types,
     * `geo::Point_t` and `geo::Vector_t`.
     * 
     * Remember the `LocalPoint_t` and `LocalVector_t` vectors from different
     * instances of `geo::OpDetGeo` have the same type but are not compatible.
     */
    
    /// Tag for vectors in the "local" GDML coordinate frame of the TPC.
    struct OpDetGeoCoordinatesTag {};
    
    /// Type of points in the local GDML TPC frame.
    using LocalPoint_t = geo::Point3DBase_t<OpDetGeoCoordinatesTag>;
    
    /// Type of displacement vectors in the local GDML TPC frame.
    using LocalVector_t = geo::Vector3DBase_t<OpDetGeoCoordinatesTag>;
    
    ///@}
    
    OpDetGeo(TGeoNode const& node, geo::TransformationMatrix&& trans);

    void   GetCenter(double* xyz, double localz=0.0) const;
    geo::Point_t const& GetCenter() const { return fCenter; }
    double RMin() const;
    double RMax() const;
    double HalfL() const;
    double HalfW() const;
    double HalfH() const;
    double Length() const { return 2.0 * HalfL(); }
    double Width() const { return 2.0 * HalfW(); }
    double Height() const { return 2.0 * HalfH(); }
    double ThetaZ() const;  ///< returns angle of detector
                            ///< with respect to z axis 
                            ///< in the Y-Z plane, in radians
    double ThetaZ(bool degrees) const; ///< returns angle of detector
                                       ///< with respect to z axis 
                                       ///< in the Y-Z plane
    //@{
    /// Get cos(angle) to normal of this detector - used for solid angle calcs
    double CosThetaFromNormal(geo::Point_t const& point) const;
    double CosThetaFromNormal(double const* xyz) const;
    //@}
    //@{
    /// Returns the distance of the specified point from detector center [cm]
    double DistanceToPoint(geo::Point_t const& point) const;
    double DistanceToPoint(double const* xyz) const;
    //@}


    /// @{
    /**
     * @name Coordinate transformation
     * 
     * Local points and displacement vectors are described by the types
     * `geo::OpDetGeo::LocalPoint_t` and `geo::OpDetGeo::LocalVector_t`,
     * respectively.
     */
    
    /// Transform point from local optical detector frame to world frame.
    void LocalToWorld(const double* opdet, double* world) const
      { fTrans.LocalToWorld(opdet, world); }
      
    /// Transform point from local optical detector frame to world frame.
    geo::Point_t toWorldCoords(LocalPoint_t const& local) const
      { return fTrans.toWorldCoords(local); }
    
    /// Transform direction vector from local to world.
    void LocalToWorldVect(const double* opdet, double* world) const
      { fTrans.LocalToWorldVect(opdet, world); }
    
    /// Transform direction vector from local to world.
    geo::Vector_t toWorldCoords(LocalVector_t const& local) const
      { return fTrans.toWorldCoords(local); }
    
    /// Transform point from world frame to local optical detector frame.
    void WorldToLocal(const double* world, double* opdet) const
      { fTrans.WorldToLocal(world, opdet); }
    
    /// Transform point from world frame to local optical detector frame.
    LocalPoint_t toLocalCoords(geo::Point_t const& world) const
      { return fTrans.toLocalCoords(world); }
    
    /// Transform direction vector from world to local.
    void WorldToLocalVect(const double* world, double* opdet) const
      { fTrans.WorldToLocalVect(world, opdet); }
    
    /// Transform direction vector from world to local.
    LocalVector_t toLocalCoords(geo::Vector_t const& world) const
      { return fTrans.toLocalCoords(world); }
    
    /// @}
    
    /// Returns the ROOT object describing the detector geometry.
    const TGeoNode*     Node() const { return fOpDetNode; }
    
    /// Returns the geometry object as `TGeoShape`.
    TGeoShape const* Shape() const { return Node()->GetVolume()->GetShape(); }
    
    /// Returns whether the detector shape is a cilynder (`TGeoTube`).
    bool isTube() const { return asTube() != nullptr; }
    
    /// Returns whether the detector shape is a bar (`TGeoBBox`).
    bool isBar() const { return (asBox() != nullptr) && !isTube(); }
    
    /**
     * @brief Prints information about this optical detector.
     * @tparam Stream type of output stream to use
     * @param out stream to send the information to
     * @param indent prepend each line with this string
     * @param verbosity amount of information printed
     * 
     * Note that the first line out the output is _not_ indented.
     * 
     * Verbosity levels
     * -----------------
     * 
     * * 0 _(default)_: only center
     * * 1: also size
     * * 2: also angle from z axis
     * 
     * The constant `MaxVerbosity` is set to the highest supported verbosity
     * level.
     */
    template <typename Stream>
    void PrintOpDetInfo
      (Stream&& out, std::string indent = "", unsigned int verbosity = 0) const;
    
    /**
     * @brief Returns a string with optical detector information
     * @see `PrintOpDetInfo()`
     * 
     * Arguments and provided information are the same as in `PrintOpDetInfo()`.
     */
    std::string OpDetInfo
      (std::string indent = "", unsigned int verbosity = 0) const;
    
    /// Maximum verbosity supported by `PrintOpDetInfo()`.
    static constexpr unsigned int MaxVerbosity = 2;
    
  private:
    using LocalTransformation_t = geo::LocalTransformationGeo
      <ROOT::Math::Transform3D, LocalPoint_t, LocalVector_t>;
    
    LocalTransformation_t fTrans; ///< Optical-detector-to-world transformation.
    const TGeoNode* fOpDetNode;  ///< Pointer to theopdet node
    geo::Point_t fCenter; ///< Stored geometric center of the optical detector.
    
    /// Returns the geometry object as `TGeoTube`, `nullptr` if not a tube.
    TGeoTube const* asTube() const
      { return dynamic_cast<TGeoTube const*>(Shape()); }
    
    /// Returns the geometry object as `TGeoBBox`, `nullptr` if not a tube.
    TGeoBBox const* asBox() const
      { return dynamic_cast<TGeoBBox const*>(Shape()); }
    
  }; // class OpDetGeo
  
} // namespace geo


//------------------------------------------------------------------------------
//--- template implementation
//---
template <typename Stream>
void geo::OpDetGeo::PrintOpDetInfo(
  Stream&& out,
  std::string indent /* = "" */,
  unsigned int verbosity /* = 0 */
) const {
  
  //----------------------------------------------------------------------------
  out << "centered at " << GetCenter() << " cm";
  
  if (verbosity-- <= 0) return; // 0
  
  //----------------------------------------------------------------------------
  if (isTube()) {
    out << ", radius: " << RMax() << " cm";
    if (RMin() != 0.0) out << " (inner: " << RMin() << " cm)";
    out << ", length: " << Length() << " cm";
  }
  else if (isBar()) {
    out << ", bar size " << Width() << " x " << Height() << " x " << Length()
      << " cm";
  }
  else out << ", shape: '" << Shape()->IsA()->GetName() << "'"; 
  
  if (verbosity-- <= 0) return; // 1
  
  //----------------------------------------------------------------------------
  out << ", theta(z): " << ThetaZ() << " rad";
  
//  if (verbosity-- <= 0) return; // 2
  
  //----------------------------------------------------------------------------
  
} // geo::OpDetGeo::PrintOpDetInfo()


#endif // LARCOREALG_GEOMETRY_OPDETGEO_H
