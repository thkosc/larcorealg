////////////////////////////////////////////////////////////////////////
/// \file  TPCGeo.h
/// \brief Encapsulate the construction of a single detector plane
///
/// \version $Id: TPCGeo.h,v 1.7 2009/12/01 21:07:51 brebel Exp $
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef GEO_TPCGEO_H
#define GEO_TPCGEO_H
#include <vector>
#include <array>

#include "TGeoVolume.h"
#include "TVector3.h"


#include "larcore/Geometry/WireGeo.h"
// TODO include PlaneGeo.h directly when issue #14365 is resolved
// #include "larcore/Geometry/Geo.h"


#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcore/Geometry/BoxBoundedGeo.h"
#include "larcore/Geometry/GeoObjectSorter.h"
#include "larcore/Geometry/LocalTransformation.h"

class TGeoNode;

namespace geo {
  class PlaneGeo;

  //......................................................................
  /// Geometry information for a single tpc
  class TPCGeo: public BoxBoundedGeo {
  public:
    
    using GeoNodePath_t = geo::WireGeo::GeoNodePath_t;
    
    static TVector3 const DirX; ///< unit vector toward positive x
    static TVector3 const DirY; ///< unit vector toward positive y
    static TVector3 const DirZ; ///< unit vector toward positive z
    
    
    // Construct a representation of a single plane of the detector
    TPCGeo(GeoNodePath_t& path, size_t depth);
    ~TPCGeo();
    
    
    /// @{
    /// @name TPC properties
    
    /// Width is associated with x coordinate [cm]
    double            ActiveHalfWidth()                         const { return fActiveHalfWidth;        }
    /// Height is associated with y coordinate [cm]
    double            ActiveHalfHeight()                        const { return fActiveHalfHeight;       }
    /// Length is associated with z coordinate [cm]
    double            ActiveLength()                            const { return fActiveLength;           }
    /// Width is associated with x coordinate [cm]
    double            HalfWidth()                               const { return fHalfWidth;              }
    /// Width is associated with x coordinate [cm]
    double            Width()                                   const { return 2.0 * HalfWidth();       }
    /// Height is associated with y coordinate [cm]
    double            HalfHeight()                              const { return fHalfHeight;             }
    /// Height is associated with y coordinate [cm]
    double            Height()                                  const { return 2.0 * HalfHeight();      }
    /// Length is associated with z coordinate [cm]
    double            Length()                                  const { return fLength;                 }
    double            ActiveMass()                              const { return fActiveVolume->Weight(); }
    const TGeoVolume* ActiveVolume()                            const { return fActiveVolume;           }
    const TGeoVolume* TotalVolume()                             const { return fTotalVolume;            }
    
    /// Returns the direction Width() is measured on
    TVector3 const&   WidthDir()                                const { return fWidthDir;               }
    
    /// Returns the direction Height() is measured on
    TVector3 const&   HeightDir()                               const { return fHeightDir;              }
    
    /// Returns the direction Length() is measured on
    TVector3 const&   LengthDir()                               const { return fLengthDir;              }
    
    DriftDirection_t  DriftDirection()                          const { return fDriftDirection;         }
    
    /// Returns the direction of the drift (vector pointing toward the planes)
    TVector3 const&   DriftDir()                                const { return fDriftDir;               }
    
    /// Drift distance is defined as the distance between the last anode plane
    /// and the opposite face of the TPC, in centimeters.
    double            DriftDistance()                           const { return ComputeDriftDistance();  }
    
    /// @}
    
    
    /// @{
    /// @name Plane access
    
    //@{
    /// Number of planes in this tpc
    unsigned int    Nplanes()                                   const { return fPlanes.size();   }
    unsigned int    NElements()                                 const { return fPlanes.size();   }
    //@}
    
    //@{
    /**
     * @brief Returns whether a plane with index iplane is present in this TPC
     * @param iplane index of plane in this TPC
     * @return whether the plane with index iplane is present in this TPC
     */
    bool HasPlane(unsigned int iplane) const { return iplane < Nplanes(); }
    bool HasElement(unsigned int iplane) const { return HasPlane(iplane); }
    //@}
    
    //@{
    /**
     * @brief Returns whether the plane in planeid is present in this TPC
     * @param planeid full plane ID
     * @return whether the plane in planeid is present in this TPC
     *
     * The cryostat and TPC numbers in planeid are ignored, as it is ignored
     * whether planeid is invalid.
     */
    bool HasPlane(geo::PlaneID const& planeid) const
      { return HasPlane(planeid.Plane); }
    bool HasElement(geo::PlaneID const& planeid) const
      { return HasPlane(planeid); }
    //@}
    
    /// Return the plane in the tpc with View_t view.
    PlaneGeo const&     Plane(geo::View_t  view)                const;
    
    
    /// Return the iplane'th plane in the TPC.
    /// @throws cet::exception (category "PlaneOutOfRange")  if no such plane
    PlaneGeo const&     Plane(unsigned int iplane)              const;
    
    //@{
    /**
     * @brief Returns the plane in planeid from this TPC
     * @param planeid full plane ID
     * @return a constant reference to the plane in planeid
     * @throws cet::exception (category "PlaneOutOfRange") if no such plane
     *
     * The cryostat and TPC numbers in planeid are ignored, as it is ignored
     * whether planeid is invalid.
     */
    const PlaneGeo&     Plane(PlaneID const& planeid)           const
      { return Plane(planeid.Plane); }
    const PlaneGeo&     GetElement(PlaneID const& planeid)            const
      { return Plane(planeid); }
    //@}
    
    /**
     * @brief Returns the plane number iplane from this TPC
     * @param iplane the number of local plane
     * @return a constant pointer to the plane, or nullptr if it does not exist
     */
    PlaneGeo const*     PlanePtr(unsigned int iplane)                const
      { return HasPlane(iplane)? fPlanes[iplane]: nullptr; }
    
    //@{
    /**
     * @brief Returns the plane in planeid from this TPC
     * @param planeid full plane ID
     * @return a constant pointer to the plane, or nullptr if it does not exist
     *
     * The cryostat and TPC numbers in planeid are ignored, as it is ignored
     * whether planeid is invalid.
     */
    PlaneGeo const*     PlanePtr(PlaneID const& planeid)                   const
      { return PlanePtr(planeid.Plane); }
    PlaneGeo const*     GetElementPtr(PlaneID const& planeid)            const
      { return PlanePtr(planeid); }
    //@}
    
    /// Returns the wire plane with the smallest surface
    geo::PlaneGeo const& SmallestPlane() const;
    
    /// Returns the first wire plane
    geo::PlaneGeo const& FirstPlane() const { return *(fPlanes[0]); }
    
    /// @brief Returns the largest number of wires among the planes in this TPC
    unsigned int MaxWires() const;
    
    /// @}
    
    /// @{
    /// @name TPC geometry properties
    
    /// Returns the center of the TPC volume in world coordinates [cm]
    TVector3 GetCenter() const;
    
    /// Returns the center of the TPC active volume in world coordinates [cm]
    TVector3 const& GetActiveVolumeCenter() const { return fActiveCenter; }
    
    /// Returns the center of the active volume face opposite to the wire planes
    /// [cm]
    TVector3 GetCathodeCenter() const;
    
    /// Returns the coordinates of the center of the specified plane [cm]
    const double*     PlaneLocation(unsigned int p)             const; 
    double            Plane0Pitch(unsigned int p)               const;
    double            PlanePitch(unsigned int p1=0,
                                 unsigned int p2=1)             const;
    double            WirePitch(unsigned int w1=0,
                                unsigned int w2=1,
                                unsigned int p=0)               const;
    
    /// Returns the identifier of this TPC
    geo::TPCID const& ID() const { return fID; }
    
    /// @}
    
#if 0
    // TODO this can be enabled once issue #14365 is resolved
    
    /// @{
    /// @name Projection on a wire plane
    /// 
    /// These methods deal with projection of points and vectors on a plane,
    /// using a geometric reference base which is not dependent on the wire
    /// direction. Technically, the objects are projected on the reference
    /// plane, that happens to be the first wire plane.
    /// In practice, the important bit is that all entities are consistently
    /// decomposed into a drift component and a remaining component. If which
    /// plane the projection happens on should matter, then geo::PlaneGeo can be
    /// used directly (but note that geo::PlaneGeo defines _two_ different
    /// frames, and the names of the frame equivalent to the one used here are
    /// different).
    /// 
    
    /// Returns the plane used for reference by projection methods
    geo::PlaneGeo const& ReferencePlane() const { return FirstPlane(); }
    
    /**
     * @brief Return the direction of reference plane width.
     * @see RefDepthDir(), DriftDir()
     * 
     * The precise definition of the vector is arbitrary, but it is defined
     * to lie on the wire plane and so that RefWidthDir(), RefDepthDir() and
     * a vector opposite to DriftDir() make a orthonormal base.
     * That base (width, depth, normal) is guaranteed to be positive defined.
     */
    TVector3 const& RefWidthDir() const { return ReferencePlane().WidthDir(); }
    
    /**
     * @brief Return the direction of reference plane depth.
     * @see RefWidthDir(), DriftDir()
     * 
     * The precise definition of the vector is arbitrary, but it is defined
     * to lie on the wire plane and so that RefWidthDir(), RefDepthDir() and
     * a vector opposite to DriftDir() make a orthonormal base.
     * That base (width, depth, normal) is guaranteed to be positive defined.
     */
    TVector3 const& RefDepthDir() const { return ReferencePlane().DepthDir(); }
    
    
    /**
     * @brief Returns the distance of the specified point from reference plane
     * @param point a point in world coordinates [cm]
     * @return the signed distance from the plane
     * 
     * The distance is defined positive if the point lies in the inner side
     * with respect to the reference plane.
     * 
     * It should match the drift distance from that plane.
     */
    double DistanceFromReferencePlane(TVector3 const& point) const
      { return ReferencePlane().DistanceFromPlane(point); }
    
    
    /**
     * @brief Decomposes a 3D point in two components
     * @param point the point to be decomposed
     * @return the two components of point, on the plane and orthogonal to it
     * 
     * The point is decomposed in:
     * 
     * 1. a component orthogonal to the reference plane, expressed as a signed
     *    real number
     * 2. a component lying on the reference plane, expressed as a 2D vector
     * 
     * The distance is from the reference plane (DistanceFromReferencePlane()).
     * The projection on the plane is obtained following the same convention
     * as PointProjection().
     */
    geo::PlaneGeo::DecomposedPoint_t DecomposePoint(TVector3 const& point) const
      { return ReferencePlane().DecomposePointWidthDepth(point); }
    
    /**
     * @brief Returns the reference point used by PointProjection()
     * 
     * The returned point is such that its decomposition results in a null
     * projection and a 0 distance from the plane.
     */
    TVector3 ProjectionReferencePoint() const
      { return ReferencePlane().GetCenter(); }
    
    /**
     * @brief Returns the projection of the specified point on the plane
     * @param point the 3D point to be projected, in world coordinates
     * @return a 2D vector representing the projection of point on the plane
     * 
     * The returned vector is a 2D vector expressing the projection of the point
     * (from world coordinates) on the reference plane.
     * The vector is expressed as @f$ ( w, d ) @f$, components following the
     * width direction (`WidthDir()`) and the depth direction (`DepthDir()`)
     * respectively. The origin point is returned by ProjectionReferencePoint().
     */
    geo::PlaneGeo::DecomposedPoint_t::Projection_t PointProjection
      (TVector3 const& point) const
      { return ReferencePlane().PointWidthDepthProjection(point); }
    
    /**
     * @brief Returns the projection of the specified vector on the plane
     * @param v the 3D vector to be projected, in world units
     * @return a 2D vector representing the projection of v on the plane
     * 
     * The returned vector is a 2D vector expressing the projection of the
     * vector (from world units) on the reference plane.
     * The vector is expressed as @f$ ( w, d ) @f$, components following the
     * width direction (`WidthDir()`) and the depth direction (`DepthDir()`)
     * respectively.
     */
    geo::PlaneGeo::DecomposedPoint_t::Projection_t VectorProjection
      (TVector3 const& point) const
      { return ReferencePlane().VectorWidthDepthProjection(point); }
    
    /**
     * @brief Returns the 3D vector from composition of projection and distance
     * @param decomp decomposed point
     * @return the 3D vector from composition of projection and distance
     * @see DecomposePoint(),
     *      ComposePoint(double, DecomposedPoint_t::Projection_t const&)
     * 
     * See `ComposePoint(double, DecomposedPoint_t::Projection_t const&)`
     * for details.
     */
    TVector3 ComposePoint(DecomposedPoint_t const& decomp) const
      { return ReferencePlane().ComposePointWidthDepth(decomp); }
    
    /**
     * @brief Returns the 3D point from composition of projection and distance
     * @param distance distance of the target point from the reference plane
     * @param proj projection of the target point on the reference plane
     * @return the 3D vector from composition of projection and distance
     * @see DecomposePoint()
     * 
     * The returned vector is the sum of two 3D vectors:
     * 
     * 1. a vector parallel to the plane normal, with norm the input distance
     * 2. a vector lying on the plane, whose projection via
     *    `PointProjection()` gives the input projection
     * 
     * Given the arbitrary definition of the projection reference, it is assumed
     * that the same convention is used as in PointProjection() and
     * DecomposePoint().
     * 
     */
    TVector3 ComposePoint
      (double distance, DecomposedPoint_t::Projection_t const& proj) const
      { return ReferencePlane().ComposePointWidthDepth(distance, proj); }
    
    /// @}
    
#endif // 0
    
    /// @{
    /// @name Coordinate transformation
    
    /// Transform point from local TPC frame to world frame
    void LocalToWorld(const double* tpc, double* world) const
      { fTrans.LocalToWorld(tpc, world); }
      
    /// Transform point from local TPC frame to world frame
    TVector3 LocalToWorld(const TVector3& local) const
      { return fTrans.LocalToWorld(local); }
    
    /// Transform direction vector from local to world
    void LocalToWorldVect(const double* tpc, double* world) const
      { fTrans.LocalToWorldVect(tpc, world); }
    
    /// Transform direction vector from local to world
    TVector3 LocalToWorldVect(const TVector3& local) const
      { return fTrans.LocalToWorldVect(local); }
    
    /// Transform point from world frame to local TPC frame
    void WorldToLocal(const double* world, double* tpc) const
      { fTrans.WorldToLocal(world, tpc); }
    
    /// Transform point from world frame to local TPC frame
    TVector3 WorldToLocal(TVector3 const& world) const
      { return fTrans.WorldToLocal(world); }
    
    /// Transform direction vector from world to local
    void WorldToLocalVect(const double* world, double* tpc) const
      { fTrans.WorldToLocalVect(world, tpc); }
    
    /// Transform direction vector from world to local
    TVector3 WorldToLocalVect(TVector3 const& world) const
      { return fTrans.WorldToLocalVect(world); }
    
    /// @}
    
    
    /**
     * @brief Returns the expected drift direction based on geometry
     * 
     * The return value is coded as follow:
     * 
     * * +1: positive x
     * * +2: positive y
     * * +3: positive z
     * * -1: negative x
     * * -2: negative y
     * * -3: negative z
     * *  0: other (or algorithm failed)
     * 
     * The current implementation is based on the assumption that electrons in
     * the middle of TPC will drift toward the wire planes, and it "never
     * fails".
     */
    short int DetectDriftDirection() const;
    
    
    
    /// Apply sorting to the PlaneGeo objects
    void              SortSubVolumes(geo::GeoObjectSorter const& sorter);
    
    /// Performs all updates after cryostat has sorted TPCs
    void UpdateAfterSorting(geo::TPCID tpcid);
    
    
    /**
     * @brief Returns whether the specified coordinate is in a range
     * @param c the coordinate
     * @param min lower boundary of the range
     * @param max upper boundary of the range
     * @param wiggle expansion factor for the range
     * @return whether the specified coordinate is in a range
     *
     * If the wiggle is larger than 1, the range is expanded by the wiggle factor.
     * If the wiggle is less than 1, the range is shrinked.
     */
    static bool CoordinateContained
      (double c, double min, double max, double wiggle = 1.)
      {
        return (c >= (min > 0? min / wiggle: min * wiggle))
          && (c <= (max < 0? max / wiggle: max * wiggle));
      } // CoordinateContained()
    
    static bool CoordinateContained
      (double c, double const* range, double wiggle = 1.)
      { return CoordinateContained(c, range[0], range[1], wiggle); }
    
    
  private:
    
    void FindPlane(GeoNodePath_t& path, size_t depth);
    void MakePlane(GeoNodePath_t& path, size_t depth);
    
  private:
    using LocalTransformation_t = geo::LocalTransformation<TGeoHMatrix>;
    
    LocalTransformation_t              fTrans;          ///< TPC-to-world transformation
    
    std::vector<PlaneGeo*>             fPlanes;         ///< List of planes in this plane
    TGeoVolume*                        fActiveVolume;   ///< Active volume of LAr, called volTPCActive in GDML file 
    TGeoVolume*                        fTotalVolume;    ///< Total volume of TPC, called volTPC in GDML file
    DriftDirection_t                   fDriftDirection; ///< Direction of the electron drift in the TPC
    std::vector<double>                fPlane0Pitch;    ///< Pitch between planes
    std::vector< std::vector<double> > fPlaneLocation;  ///< xyz locations of planes in the TPC
    TVector3                           fActiveCenter;   ///< center of the active volume, in world coordinates [cm]

    double                             fActiveHalfWidth;  ///< half width of active volume
    double                             fActiveHalfHeight; ///< half height of active volume
    double                             fActiveLength;     ///< length of active volume
    double                             fHalfWidth;        ///< half width of total volume
    double                             fHalfHeight;       ///< half height of total volume
    double                             fLength;           ///< length of total volume

    TVector3 fWidthDir; ///< direction width refers to
    TVector3 fHeightDir; ///< direction height refers to
    TVector3 fLengthDir; ///< direction length refers to
    TVector3 fDriftDir; ///< direction electrons drift along
    
    geo::TPCID                         fID;             ///< ID of this TPC
    
    /// Index of the plane for each view (InvalidID if none)
    std::vector<geo::PlaneID::PlaneID_t> fViewToPlaneNumber;
  
    /// Recomputes the drift direction; needs planes to have been initialised
    void ResetDriftDirection();
    
    /// Computes the distance between the cathode and the last wire plane
    /// (last respect to the sorting order)
    double ComputeDriftDistance() const;
    
    /// Recomputes the TPC boundary
    void InitTPCBoundaries();
  
  };
}

#endif
////////////////////////////////////////////////////////////////////////
