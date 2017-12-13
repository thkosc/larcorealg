////////////////////////////////////////////////////////////////////////
/// \file  larcorealg/Geometry/PlaneGeo.h
/// \brief Encapsulate the construction of a single detector plane
///
/// \author  brebel@fnal.gov
////////////////////////////////////////////////////////////////////////
#ifndef LARCOREALG_GEOMETRY_PLANEGEO_H
#define LARCOREALG_GEOMETRY_PLANEGEO_H

// LArSoft libraries
#include "larcorealg/CoreUtils/DereferenceIterator.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/Geometry/SimpleGeo.h"
#include "larcorealg/Geometry/LocalTransformationGeo.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/Decomposer.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/geo_vectors_utils.h" // geo::vect
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// ROOT libraries
#include "Math/GenVector/Cartesian2D.h"
#include "Math/GenVector/PositionVector2D.h"
#include "Math/GenVector/DisplacementVector2D.h"

// ROOT libraries
#include "TGeoMatrix.h" // TGeoHMatrix

// C/C++ standard libraries
#include <cmath> // std::atan2()
#include <vector>
#include <string>


class TGeoNode;
class TVector3;

namespace geo {

  namespace details {
    class ActiveAreaCalculator;
  } // namespace details
  
  //......................................................................
  
  /**
   * @brief Geometry information for a single wire plane.
   * 
   * The plane is represented in the geometry by a solid which contains wires.
   * Currently, only box solids are well supported.
   * The box which is representation of the plane has some thickness, and it
   * should not be assumed that the wires are in the median section of it,
   * that is, the center of the box may not lie on the plane defined by the
   * wires.
   * 
   * The plane defines two local reference frames.
   * The first, depending on wire directions and therefore called "wire base",
   * is defined by the normal to the plane (pointing toward the center of the
   * TPC), the direction of the wires, and the direction that the wires measure.
   * This is a positive orthogonal base.
   * Note that for this base to be correctly defined, the Geometry service has
   * to provide external information (for example, where the center of the
   * TPC is).
   * 
   * The second, depending only on the shape of the plane and called "frame
   * base", is defined by the normal (the same as for the previous one), and two
   * orthogonal axes, "width" and "depth", aligned with the sides of the plane.
   * If the plane has not the shape of a box, this reference frame is not
   * available. This coordinate system is also positive defined.
   * These components are all measured in centimeters.
   * 
   */
  // Note: SignalType() and SetSignalType() have been removed.
  //       Use `geo::GeometryCore::SignalType` instead.
  //       (see LArSoft issue #14365 at https://cdcvs.fnal.gov/redmine/issues/14365 )
  class PlaneGeo {
    
    using DefaultVector_t = TVector3; // ... not for long
    using DefaultPoint_t = TVector3; // ... not for long
    
  public:
    
    using GeoNodePath_t = std::vector<TGeoNode const*>;
    
    /// @{
    /**
     * @name Types for geometry-local reference vectors.
     * 
     * These types represents points and displacement vectors in the reference
     * frame defined in the plane geometry box from the GDML geometry
     * description.
     * 
     * No alias is explicitly defined for the LArSoft global vector types,
     * `geo::Point_t` and `geo::Vector_t`.
     * 
     * Remember the `LocalPoint_t` and `LocalVector_t` vectors from different
     * instances of `geo::PlaneGeo` have the same type but are not compatible.
     */
    
    /// Tag for vectors in the "local" GDML coordinate frame of the plane.
    struct PlaneGeoCoordinatesTag {};
    
    /// Type of points in the local GDML wire plane frame.
    using LocalPoint_t = geo::Point3DBase_t<PlaneGeoCoordinatesTag>;
    
    /// Type of displacement vectors in the local GDML wire plane frame.
    using LocalVector_t = geo::Vector3DBase_t<PlaneGeoCoordinatesTag>;
    
    ///@}
    
    /// @{
    /// @name Types for vectors in the wire coordinate frame.
    
    /// Tag for wire base vectors.
    struct WireCoordinateReferenceTag {};
    
    /// Type for projections in the wire base representation.
    using WireCoordProjection_t = ROOT::Math::DisplacementVector2D
      <ROOT::Math::Cartesian2D<double>, WireCoordinateReferenceTag>;
    
    /// Type used for plane decompositions on wire base.
    using WireDecomposer_t = geo::Decomposer
      <geo::Vector_t, geo::Point_t, WireCoordProjection_t>;
    
    /// Type describing a 3D point or vector decomposed on a plane on wire base.
    using WireDecomposedVector_t = WireDecomposer_t::DecomposedVector_t;
    
    /// @}
    
    
    /// @{
    /// @name Types for vectors in the width/depth coordinate frame.
    
    /// Tag for plane frame base vectors.
    struct WidthDepthReferenceTag {};
    
    /// Type for projections in the plane frame base representation.
    /// @fixme the following should be a PositionVector2D
    using WidthDepthProjection_t = ROOT::Math::DisplacementVector2D
      <ROOT::Math::Cartesian2D<double>, WidthDepthReferenceTag>;
    
    /// Type for vector projections in the plane frame base representation.
    using WidthDepthDisplacement_t = ROOT::Math::DisplacementVector2D
      <ROOT::Math::Cartesian2D<double>, WidthDepthReferenceTag>;
    
    /// Type used for plane decompositions on plane frame (width/depth).
    using WidthDepthDecomposer_t = geo::Decomposer
      <geo::Vector_t, geo::Point_t, WidthDepthProjection_t>;
    
    /// Type describing a 3D point or vector decomposed on a plane
    /// with plane frame base (width and depth).
    using WDDecomposedVector_t = WidthDepthDecomposer_t::DecomposedVector_t;
    
    /// @}
    
    
    /// Type for description of rectangles.
    using Rect = lar::util::simple_geo::Rectangle<double>;
    
    
    /// Construct a representation of a single plane of the detector
    PlaneGeo(GeoNodePath_t& path, size_t depth);
    
    
    /// @{
    /// @name Plane properties
    
    /// Which coordinate does this plane measure
    View_t View()                                             const { return fView;          }
    
    /// What is the orientation of the plane
    Orient_t Orientation()                                    const { return fOrientation;   }

    /// Angle of the wires from positive z axis; @f$ \theta_{z} \in [ 0, \pi ]@f$.
    double ThetaZ()                                           const;
    
    /// Angle from positive z axis of the wire coordinate axis, in radians
    double PhiZ()                                             const
      { return std::atan2(fSinPhiZ, fCosPhiZ); }
    
    /// Sine of PhiZ()
    double SinPhiZ()                                          const { return fSinPhiZ; }
    
    /// Cosine of PhiZ()
    double CosPhiZ()                                          const { return fCosPhiZ; }
    
    /// Returns the identifier of this plane
    geo::PlaneID const& ID() const { return fID; }
    
    /// @}
    
    /// @{
    /// @name Plane size and coordinates
    
    //@{
    /**
     * @brief Return the direction of plane width.
     * @tparam Vector the type of vector to return (current default: `TVector3`)
     * 
     * The precise definition of the sides is arbitrary, but they are defined
     * to lie on the wire plane and so that WidthDir(), DepthDir() and
     * GetNormalDirection() make a orthonormal base.
     * That base (width, depth, normal) is guaranteed to be positive defined.
     */
    template <typename Vector = DefaultVector_t>
    Vector WidthDir() const { return geo::vect::convertTo<Vector>(fDecompFrame.MainDir()); }
    //@}
    
    /**
     * @brief Return the direction of plane depth.
     * @tparam Vector the type of vector to return (current default: `TVector3`)
     * 
     * The precise definition of the sides is arbitrary, but they are defined
     * to lie on the wire plane and so that WidthDir(), DepthDir() and
     * GetNormalDirection() make a orthonormal base.
     * That base (width, depth, normal) is guaranteed to be positive defined.
     */
    template <typename Vector = DefaultVector_t>
    Vector DepthDir() const { return geo::vect::convertTo<Vector>(fDecompFrame.SecondaryDir()); }
    
    /**
     * @brief Return the width of the plane.
     * @see Depth(), WidthDir(), DepthDir()
     * 
     * The precise definition is arbitrary (see `WidthDir()`).
     */
    double Width() const { return fFrameSize.Width(); }
    
    /**
     * @brief Return the depth of the plane.
     * @see Width(), WidthDir(), DepthDir()
     * 
     * The precise definition is arbitrary (see `DepthDir()`).
     */
    double Depth() const { return fFrameSize.Depth(); }

    
    /// Returns the world coordinates of the box containing the plane.
    /// @see GetBoxCenter()
    geo::BoxBoundedGeo BoundingBox() const;
    
    /// @}
    
    
    /// @{
    /// @name Wire access
    
    //@{
    /// Number of wires in this plane
    unsigned int Nwires()                                     const { return fWire.size();   }
    unsigned int NElements()                                  const { return Nwires();       }
    //@}
    
    //@{
    /**
     * @brief Returns whether a wire with index iwire is present in this plane.
     * @param iwire index of wire in this plane
     * @return whether the wire with index iwire is present in this plane
     */
    bool HasWire(unsigned int iwire) const { return iwire < Nwires(); }
    bool HasElement(unsigned int iwire) const { return HasWire(iwire); }
    //@}
    
    //@{
    /**
     * @brief Returns whether the wire in wireid is present in this plane.
     * @param wireid full wire ID
     * @return whether the wire in wireid is present in this plane
     *
     * The cryostat, TPC and plane numbers in wireid are ignored, as it is
     * ignored whether wireid is invalid.
     */
    bool HasWire(geo::WireID const& wireid) const
      { return HasWire(wireid.Wire); }
    bool HasElement(geo::WireID const& wireid) const
      { return HasWire(wireid); }
    //@}
    
    /// Return the iwire'th wire in the plane.
    /// @throws cet::exception (category "WireOutOfRange") if no such wire
    /// @note In the past, no check was performed.
    WireGeo const& Wire(unsigned int iwire) const;
    
    //@{
    /**
     * @brief Returns the wire in wireid from this plane.
     * @param wireid full wire ID
     * @return a constant reference to the wire in wireid
     * @throws cet::exception (category "WireOutOfRange") if no such wire
     *
     * The cryostat, TPC and plane numbers in wireid are ignored, as it is
     * ignored whether wireid is invalid.
     */
    WireGeo const& Wire(WireID const& wireid) const
      { return Wire(wireid.Wire); }
    WireGeo const& GetElement(WireID const& wireid) const
      { return Wire(wireid); }
    //@}
    
    /**
     * @brief Returns the wire number iwire from this plane.
     * @param iwire the number of local wire
     * @return a constant pointer to the wire, or nullptr if it does not exist
     */
    WireGeo const* WirePtr(unsigned int iwire) const
      { return HasWire(iwire)? &(fWire[iwire]): nullptr; }
    
    //@{
    /**
     * @brief Returns the wire in wireid from this plane.
     * @param wireid full wire ID
     * @return a constant pointer to the wire, or nullptr if it does not exist
     *
     * The cryostat, TPC and plane numbers in wireid are ignored, as it is
     * ignored whether wireid is invalid.
     */
    WireGeo const* WirePtr(WireID const& wireid) const
      { return WirePtr(wireid.Wire); }
    WireGeo const* GetElementPtr(WireID const& wireid) const
      { return WirePtr(wireid); }
    //@}
    
    
    /// Return the first wire in the plane.
    const WireGeo& FirstWire()                                const { return Wire(0);        }
    
    /// Return the middle wire in the plane.
    const WireGeo& MiddleWire()                               const { return Wire(Nwires()/2); }
    
    /// Return the last wire in the plane.
    const WireGeo& LastWire()                                 const { return Wire(Nwires()-1); }
    
    /**
     * @brief Allows range-for iteration on all wires in this plane.
     * @return an object suitable for range-for iteration on all wires
     * 
     * This example uses geometry to iterate on all planes in an outer loop,
     * and then iterates on all wires in the plane in the inner loop:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * 
     * auto geom = lar::providerFrom<geo::Geometry>();
     * for (geo::PlaneGeo const& plane: geom->IteratePlanes()) {
     *   
     *   // collect plane information
     *   
     *   for (geo::WireGeo const& wire: plane.IterateWires()) {
     *     
     *     // do something with each single wire
     *     
     *   }
     * } // for planes
     * 
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * (note that all data types here can be replaced with `auto`).
     * 
     */
    auto IterateWires() const
      { return fWire; }
    
    /// @}
    
    
    /// @{
    /// @name Plane geometry properties
    
    /// Return the wire pitch (in centimeters). It is assumed constant.
    double WirePitch() const { return fWirePitch; }
    
    /**
     * @brief Returns whether the higher z wires have higher wire ID.
     * @return whether the higher z wires have higher wire ID
     * @see GetIncreasingWireDirection()
     * 
     * This method is related to GetIncreasingWireDirection()
     * (it might be expressed as "GetIncreasingWireDirection()[2] > 0"),
     * but it is implemented in a faster and independent way.
     */
    bool WireIDincreasesWithZ() const;
    
    /**
     * @brief Returns the direction normal to the plane.
     * @tparam Vector the type of vector to return (current default: `TVector3`)
     * @return a TVector3 versor with a direction normal to the plane
     *
     * The versor is orthogonal to the plane.
     * The direction is defined so that the semi-space pointed to contains
     * the TPC center.
     * 
     * @note Each decomposition base (wire-based and frame-based) has its own
     *       normal, defined solely from its two decomposition plane axes.
     *       The wire-based frame is nevertheless required to have a normal
     *       matching this one, while the frame-based normal might happen to be
     *       in the opposite direction depending on the original geometry
     *       description.
     */
    template <typename Vector = DefaultVector_t>
    Vector GetNormalDirection() const { return geo::vect::convertTo<Vector>(fNormal); }

    /**
     * @brief Returns the direction of increasing wires.
     * @tparam Vector the type of vector to return (current default: `TVector3`)
     * @return a TVector3 versor with the direction of increasing wires
     *
     * The versor is orthogonal to the wires (assumed parallel),
     * lies on the plane and its direction goes toward increasing wire IDs.
     */
    template <typename Vector = DefaultVector_t>
    Vector GetIncreasingWireDirection() const
      { return geo::vect::convertTo<Vector>(fDecompWire.SecondaryDir()); }
    
    
    /**
     * @brief Returns the centre of the wire plane in world coordinates [cm]
     * @see GetBoxCenter()
     * 
     * The center of the plane is defined so that it has width and depth
     * coordinates in the middle of the plane box (that is, the geometrical
     * representation of the plane in the geometry description), and the other
     * coordinate set at drift distance 0.
     * 
     * Note that this does not necessarily match the center of the box, if the
     * geometry does not place the wires, which define the drift distance, in
     * the plane in the middle of the box.
     */
    template <typename Point = DefaultPoint_t>
    Point GetCenter() const
      { return geo::vect::convertTo<Point>(fCenter); }
    
    /**
     * @brief Returns the centre of the box representing the plane.
     * @tparam Point type of point to be returned (current default: `TVector3`)
     * @return the centre of the box, in world coordinates [cm]
     * @see GetCenter()
     * 
     * This is the centre of the box representing the plane in the geometry
     * description, in world coordinates.
     * This is rarely of any use, as most of the times `GetCenter()` delivers
     * the proper information, e.g. for simulation and reconstruction.
     */
    template <typename Point = DefaultPoint_t>
    Point GetBoxCenter() const
      { return geo::vect::convertTo<Point>(toWorldCoords(LocalPoint_t{ 0.0, 0.0, 0.0 })); }
    
    /**
     * @brief Returns the direction of the wires.
     * @tparam Vector the type of vector to return (current default: `TVector3`)
     * @return a unit vector following the direction of the wires
     * 
     * All wires in the plane are assumed parallel.
     */
    template <typename Vector = DefaultVector_t>
    Vector GetWireDirection() const { return geo::vect::convertTo<Vector>(fDecompWire.MainDir()); }
    
    
    //@{
    /**
     * @brief Returns the ID of wire closest to the specified position.
     * @param pos world coordinates of the point [cm]
     * @return the ID of the wire closest to the projection of pos on the plane
     * @throw InvalidWireError (category: `"Geometry"`) if out of range
     *
     * The position is projected on the wire plane, and the ID of the nearest
     * wire to the projected point is returned.
     * 
     * If the wire does not exist, an exception is thrown that contains both the
     * wire that would be the closest one (`badWireID()`), and also the wire
     * that is actually the closest one (`betterWireID()`). When this happens,
     * the specified position was outside the wire plane.
     * 
     * Note that the caller should check for containment: this function may or
     * may not report the position being outside the plane, depending on where
     * it is. In the current implementation, the wires are considered infinitely
     * long, and if the position projection is closer than half the wire pitch
     * from any of these extrapolated wires, the method will not report error.
     * 
     */
    geo::WireID NearestWireID(geo::Point_t const& pos) const;
    geo::WireID NearestWireID(TVector3 const& pos) const
      { return NearestWireID(geo::vect::toPoint(pos)); }
    //@}
    
    
    //@{
    /**
     * @brief Returns the distance of the specified point from the wire plane.
     * @param point a point in world coordinates [cm]
     * @return the signed distance from the wire plane
     * 
     * The distance is defined positive if the point lies in the side the normal
     * vector (GetNormalDirection()) points to.
     * 
     * The distance is defined from the geometric plane where the wires lie, and
     * it may not match the distance from the center of the geometry box
     * representing the plane.
     * It should always match the drift distance from this wire plane, and the
     * result of `DriftPoint(point, DistanceFromPlane(point))` will bring the
     * point to the plane.
     */
    double DistanceFromPlane(geo::Point_t const& point) const
      { return fDecompWire.PointNormalComponent(point); }
    double DistanceFromPlane(TVector3 const& point) const
      { return DistanceFromPlane(geo::vect::toPoint(point)); }
    //@}
    
    
    //@{
    /**
     * @brief Shifts the position of an electron drifted by a distance.
     * @param position _(modified)_ the position of the electron
     * @param drift drift distance to shift the electron by [cm]
     * 
     * This is a pure geometry computation: the position is shifted by the drift
     * distance in the direction opposite to the normal to the plane (as
     * returned by `GetNormalDirection()`), no matter where the position is
     * relative to the plane.
     * The wording about "electron position" is just meant to remind that the
     * drift shift is taken with opposite sign: since the point is assumed to be
     * an electron, a positive drift normally moves its position toward the wire
     * plane.
     */
    void DriftPoint(geo::Point_t& position, double distance) const
      { position -= distance * GetNormalDirection<geo::Vector_t>(); }
    void DriftPoint(TVector3& position, double distance) const
      { position -= distance * GetNormalDirection(); }
    //@}
    
    //@{
    /**
     * @brief Shifts the position along drift direction to fall on the plane.
     * @param position _(modified)_ the position to be shifted
     * 
     * This is a pure geometry computation: the position is shifted by the drift
     * distance in the direction opposite to the normal to the plane (as
     * returned by `GetNormalDirection()`), no matter where the position is
     * relative to the plane.
     */
    void DriftPoint(geo::Point_t& position) const
      { DriftPoint(position, DistanceFromPlane(position)); }
    void DriftPoint(TVector3& position) const
      { DriftPoint(position, DistanceFromPlane(position)); }
    //@}
    
    
    /**
     * @brief Returns an area covered by the wires in the plane.
     * 
     * The returned value is conceptually akin of a projection of `Coverage()`
     * volume. Yet, the precise definition of the area is not specified,
     * therefore this area should not be uses for physics.
     * 
     * The current implementation is documented in
     * `details::ActiveAreaCalculator`.
     * 
     */
    Rect const& ActiveArea() const { return fActiveArea; }
    
    /// Returns a volume including all the wires in the plane.
    lar::util::simple_geo::Volume<> Coverage() const;
    
    /**
     * @brief Prints information about this plane.
     * @tparam Stream type of output stream to use
     * @param out stream to send the information to
     * @param indent prepend each line with this string
     * @param verbosity amount of information printed
     * 
     * Information on single wires is not printed.
     * Note that the first line out the output is _not_ indented.
     * 
     * Verbosity levels
     * -----------------
     * 
     * * 0: only plane ID
     * * 1 _(default)_: also center and wire angle
     * * 2: also information about wires
     * * 3: also information about normal and increasing coordinate direction
     * * 4: also information about wire direction, width and depth
     * * 5: also coverage
     * * 6: also bounding box
     * 
     */
    template <typename Stream>
    void PrintPlaneInfo
      (Stream&& out, std::string indent = "", unsigned int verbosity = 1) const;
    
    /// Maximum value for print verbosity.
    static constexpr unsigned int MaxVerbosity = 6;
    
    /// @}
    
    
    /// @{
    /// @name Projections on wire length/wire coordinate direction base
    /// 
    /// These methods deal with projection of points and vectors on the plane,
    /// using a geometric reference base which is dependent on the wire
    /// direction. This is useful for plane reconstruction.
    /// 
    
    //@{
    /**
     * @brief Returns the coordinate of point on the plane respect to a wire.
     * @param point world coordinate of the point to get the coordinate of [cm]
     * @param refWire reference wire
     * @return the coordinate of the point [cm]
     * @see WireCoordinate()
     * 
     * The method returns the coordinate of the point in the direction measured
     * by the wires on this plane starting from the specified reference wire,
     * in world units (that is, centimeters).
     *  
     * The point does not need to be on the plane, and the projection of the
     * point to the plane is considered.
     * The reference wire, instead, must belong to this plane. This assumption
     * is not checked, and if violated the results are undefined (in the current
     * implementation, they are just wrong).
     */
    double PlaneCoordinateFrom
      (geo::Point_t const& point, geo::WireGeo const& refWire) const
      {
        return
          fDecompWire.VectorSecondaryComponent(point - geo::vect::toPoint(refWire.GetCenter()));
      }
    double PlaneCoordinateFrom
      (TVector3 const& point, geo::WireGeo const& refWire) const
      { return PlaneCoordinateFrom(geo::vect::toPoint(point), refWire); }
    //@}
    
    
    //@{
    /**
     * @brief Returns the coordinate of the point on the plane.
     * @param point world coordinate of the point to get the coordinate of [cm]
     * @return the coordinate of the point [cm]
     * @see PlaneCoordinateFrom(TVector3 const&, geo::Wire const&)
     * 
     * The method returns the coordinate of the point in the direction measured
     * by the wires on this plane starting on the first wire, in world units
     * (that is, centimeters). A point on the first wire will have coordinate
     * 0.0, one on the next wire will have coordinate equal to a single wire
     * pitch, etc.
     *  
     * The point does not need to be on the plane, and the projection of the
     * point to the plane is considered.
     */
    double PlaneCoordinate(geo::Point_t const& point) const
      { return fDecompWire.PointSecondaryComponent(point); }
    double PlaneCoordinate(TVector3 const& point) const
      { return PlaneCoordinate(geo::vect::toPoint(point)); }
    //@}
    
    /**
     * @brief Returns the coordinate of the point on the plane, in wire units.
     * @param point world coordinate of the point to get the coordinate of
     * @return the coordinate of the point, in wire pitch units
     * @see CoordinateFrom(TVector3 const&, geo::Wire const&)
     * 
     * The method returns the coordinate of the point in the direction measured
     * by the wires on this plane starting on the first wire, in wire units
     * (that is, wire pitches). A point on the first wire will have coordinate
     * 0.0, one on the next wire will have coordinate 1.0, etc.
     *  
     * The point does not need to be on the plane, and the projection of the
     * point to the plane is considered.
     */
    template <typename Point = DefaultPoint_t>
    double WireCoordinate(Point const& point) const
      { return PlaneCoordinate(point) / WirePitch(); }
    
    
    //@{
    /**
     * @brief Decomposes a 3D point in two components.
     * @param point the point to be decomposed
     * @return the two components of point, on the plane and orthogonal to it
     * 
     * The point is decomposed in:
     * 
     * 1. a component orthogonal to the plane, expressed as a signed real number
     * 2. a component lying on the plane, expressed as a 2D vector
     * 
     * The distance is obtained as by DistanceFromPlane().
     * The projection on the plane is obtained following the same convention
     * as PointProjection().
     */
    WireDecomposedVector_t DecomposePoint(geo::Point_t const& point) const
      { return fDecompWire.DecomposePoint(point); }
    WireDecomposedVector_t DecomposePoint(TVector3 const& point) const
      { return DecomposePoint(geo::vect::toPoint(point)); }
      
    //@}
    
    /**
     * @brief Returns the reference point used by `PointProjection()`.
     * @tparam Point the type of point to return (current default: `TVector3`)
     * 
     * The returned point is such that its decomposition results in a null
     * projection and a 0 distance from the plane.
     */
    template <typename Point = DefaultPoint_t>
    Point ProjectionReferencePoint() const
      { return geo::vect::convertTo<Point>(fDecompWire.ReferencePoint()); }
    
    //@{
    /**
     * @brief Returns the projection of the specified point on the plane.
     * @param point the 3D point to be projected, in world coordinates
     * @return a 2D vector representing the projection of point on the plane
     * 
     * The returned vector is a 2D vector expressing the projection of the point
     * (from world coordinates) on the wire plane.
     * The vector is expressed as @f$ ( \ell, w ) @f$. The component
     * @f$ \ell @f$ is measured on the direction of the first wire (see
     * `WireGeo::Direction()`), using its center (see `WireGeo::GetCenter()`) as
     * reference point. The component @f$ w @f$ is defined on the wire
     * coordinate direction (see `GetIncreasingWireDirection()`), relative to
     * the first wire, as it is returned by PlaneCoordinate(). 
     * 
     * The reference point is also returned by ProjectionReferencePoint().
     */
    WireCoordProjection_t Projection(geo::Point_t const& point) const
      { return fDecompWire.ProjectPointOnPlane(point); }
    WireCoordProjection_t PointProjection(geo::Point_t const& point) const
      { return Projection(point); }
    WireCoordProjection_t PointProjection(TVector3 const& point) const
      { return PointProjection(geo::vect::toPoint(point)); }
    //@}
    
    //@{
    /**
     * @brief Returns the projection of the specified vector on the plane.
     * @param v the 3D vector to be projected, in world units
     * @return a 2D vector representing the projection of v on the plane
     * 
     * The returned vector is a 2D vector expressing the projection of the
     * vector (from world units) on the wire plane.
     * The vector is expressed as @f$ ( \ell, w ) @f$. The component
     * @f$ \ell @f$ is measured on the direction of the first wire (see
     * `WireGeo::Direction()`). The component @f$ w @f$ is defined on the wire
     * coordinate direction (see `GetIncreasingWireDirection()`). 
     */
    WireCoordProjection_t Projection(geo::Vector_t const& v) const
      { return fDecompWire.ProjectVectorOnPlane(v); }
    WireCoordProjection_t VectorProjection(geo::Vector_t const& v) const
      { return Projection(v); }
    WireCoordProjection_t VectorProjection(TVector3 const& v) const
      { return VectorProjection(geo::vect::toVector(v)); }
    //@}
    
    /**
     * @brief Returns the 3D vector from composition of projection and distance.
     * @tparam Point the type of point to return (current default: `TVector3`)
     * @param decomp decomposed point
     * @return the 3D vector from composition of projection and distance
     * @see DecomposePoint(), ComposePoint(double, WireCoordProjection_t const&)
     * 
     * See `ComposePoint(double, WireCoordProjection_t const&)` for details.
     */
    template <typename Point = DefaultPoint_t>
    Point ComposePoint(WireDecomposedVector_t const& decomp) const
      { return geo::vect::convertTo<Point>(fDecompWire.ComposePoint(decomp)); }
    
    /**
     * @brief Returns the 3D point from composition of projection and distance.
     * @tparam Point the type of point to return (current default: `TVector3`)
     * @param distance distance of the target point from the wire plane
     * @param proj projection of the target point on the wire plane
     * @return the 3D vector from composition of projection and distance
     * @see DecomposePoint()
     * 
     * The returned point is the reference point of the frame system (that is,
     * the plane center), translated by two 3D vectors:
     * 
     * 1. a vector parallel to the plane normal, with norm the input distance
     * 2. a vector lying on the plane, whose projection via `PointProjection()`
     *    gives the input projection
     * 
     * The choice of the projection reference point embodies the same convention
     * used in `PointProjection()` and `DecomposePoint()`.
     * In fact, the strict definition of the result of this method is a 3D point
     * whose decomposition on the plane frame base matches the method arguments.
     * 
     * Note that currently no equivalent facility is available to compose
     * vectors instead of points, that is, entities ignoring the reference
     * point.
     */
    template <typename Point = DefaultPoint_t>
    Point ComposePoint
      (double distance, WireCoordProjection_t const& proj) const
      { return geo::vect::convertTo<Point>(fDecompWire.ComposePoint(distance, proj)); }
    
    
    /// @}
    
    
    /// @{
    /// @name Projection on width/depth plane
    /// 
    /// These methods deal with projection of points and vectors on the plane,
    /// using a geometric reference base which is not dependent on the wire
    /// direction. This is more useful when comparing with the TPC or other
    /// planes.
    /// 
    
    //@{
    /**
     * @brief Decomposes a 3D point in two components.
     * @param point the point to be decomposed
     * @return the two components of point, on the plane and orthogonal to it
     * 
     * The point is decomposed in:
     * 
     * 1. a component orthogonal to the plane, expressed as a signed real number
     * 2. a component lying on the plane, expressed as a 2D vector
     * 
     * The distance is obtained as by DistanceFromPlane().
     * The projection on the plane is obtained following the same convention
     * as PointWidthDepthProjection().
     */
    WDDecomposedVector_t DecomposePointWidthDepth(geo::Point_t const& point) const
      { return fDecompFrame.DecomposePoint(point); }
    WDDecomposedVector_t DecomposePointWidthDepth(TVector3 const& point) const
      { return DecomposePointWidthDepth(geo::vect::toPoint(point)); }
    //@}
    
    //@{
    /**
     * @brief Returns the projection of the specified point on the plane.
     * @param point the 3D point to be projected, in world coordinates
     * @return a 2D vector representing the projection of point on the plane
     * 
     * The returned vector is a 2D vector expressing the projection of the point
     * (from world coordinates) on the wire plane.
     * The vector is expressed as @f$ ( w, d ) @f$, components following the
     * width direction (`WidthDir()`) and the depth direction (`DepthDir()`)
     * respectively. The origin point is the center of the plane.
     */
    WidthDepthProjection_t PointWidthDepthProjection
      (geo::Point_t const& point) const
      { return fDecompFrame.ProjectPointOnPlane(point); }
    WidthDepthProjection_t PointWidthDepthProjection
      (TVector3 const& point) const
      { return PointWidthDepthProjection(geo::vect::toPoint(point)); }
    //@}
    
    //@{
    /**
     * @brief Returns the projection of the specified vector on the plane.
     * @param v the 3D vector to be projected, in world units
     * @return a 2D vector representing the projection of v on the plane
     * 
     * The returned vector is a 2D vector expressing the projection of the
     * vector (from world units) on the wire plane.
     * The vector is expressed as @f$ ( w, d ) @f$, components following the
     * width direction (`WidthDir()`) and the depth direction (`DepthDir()`)
     * respectively.
     */
    WidthDepthProjection_t VectorWidthDepthProjection
      (geo::Vector_t const& v) const
      { return fDecompFrame.ProjectVectorOnPlane(v); }
    WidthDepthProjection_t VectorWidthDepthProjection
      (TVector3 const& v) const
      { return VectorWidthDepthProjection(geo::vect::toVector(v)); }
    //@}
    
    //@{
    /**
     * @brief Returns if the projection of specified point is within the plane.
     * @param point world coordinate of the point to test [cm]
     * @return whether the projection of specified point is within the plane
     * @see PointWidthDepthProjection(), Width(), Height()
     * 
     * The method extracts the projection of the specified point on the plane,
     * as in `PointWidthDepthProjection()`, and then verifies that the
     * projection falls within the wire plane area, as defined by the dimensions
     * from the geometry description.
     */
    bool isProjectionOnPlane(geo::Point_t const& point) const;
    bool isProjectionOnPlane(TVector3 const& point) const
      { return isProjectionOnPlane(geo::vect::toPoint(point)); }
    //@}
    
    /**
     * @brief Returns a projection vector that, added to the argument, gives a
     *        projection inside (or at the border of) the plane.
     * @param proj starting projection
     * @param wMargin the point is brought this amount _inside_ the target area
     * @param dMargin the point is brought this amount _inside_ the target area
     * @return a projection displacement
     * @see `DeltaFromActivePlane()`
     * 
     * The returned projection vector is guaranteed, when added to `proj`, to
     * yield a projection on or within the border of the plane (the "target
     * area"), as defined by the GDML geometry.
     * 
     * The target plane area is reduced on each side by the specified margins.
     * If for example `wMargin` is `1.0`, the area lower border on the width
     * direction will be increased by 1 cm, and the upper border will be
     * decreased by 1 cm effectively making the area 2 cm narrowed on the width
     * direction.
     * The same independently applies to the depth direction with `dMargin`.
     * The main purpose of the margins is to accommodate for rounding errors.
     * A version of this method with default margins of 0 is also available.
     * 
     * If the projection is already on the target area, the returned
     * displacement is null.
     */
    WidthDepthProjection_t DeltaFromPlane
      (WidthDepthProjection_t const& proj, double wMargin, double dMargin)
       const;
    
    /**
     * @brief Returns a projection vector that, added to the argument, gives a
     *        projection inside (or at the border of) the area of plane.
     * @param proj starting projection
     * @param margin the point is brought this amount _inside_ the plane area
     *        _(default: 0)_
     * @return a projection displacement
     * @see `DeltaFromPlane(WidthDepthProjection_t const&, double, double)`
     * 
     * This is the implementation with default values for margins of
     * `DeltaFromPlane()`.
     * The depth and width margins are the same, and 0 by default.
     */
    WidthDepthProjection_t DeltaFromPlane
      (WidthDepthProjection_t const& proj, double margin = 0.0) const
      { return DeltaFromPlane(proj, margin, margin); }
    
    /**
     * @brief Returns a projection vector that, added to the argument, gives a
     *        projection inside (or at the border of) the active area of plane.
     * @param proj starting projection
     * @param wMargin the point is brought this amount _inside_ the active area
     * @param dMargin the point is brought this amount _inside_ the active area
     * @return a projection displacement
     * @see `DeltaFromPlane()`
     * 
     * The "active" area of the plane is the rectangular area which includes all
     * the wires. The area is obtained as the smallest rectangle including
     * the projection of both ends of all wires in the plane, less half a pitch.
     * This defines a "fiducial" area away from the borders of the plane.
     * The projection is in the frame reference (`PointWidthDepthProjection()`).
     * The area is reduced on each side by the specified margins. If for example
     * `wMargin` is `1.0`, the active area lower border on the width direction
     * will be increased by 1 cm, and the upper border will be decreased by 1 cm
     * effectively making the active area 2 cm narrowed on the width direction.
     * The same independently applies to the depth direction with `dMargin`.
     * The main purpose of the margins is to accommodate for rounding errors.
     * A version of this method with default margins of 0 is also available.
     * 
     * If the projection is already on the active area of the plane, the
     * returned displacement is null.
     * Otherwise, the displacement, added to proj, will bring it on the active
     * plane area (in fact, on its border).
     */
    WidthDepthProjection_t DeltaFromActivePlane
      (WidthDepthProjection_t const& proj, double wMargin, double dMargin)
      const;
    
    /**
     * @brief Returns a projection vector that, added to the argument, gives a
     *        projection inside (or at the border of) the active area of plane.
     * @param proj starting projection
     * @param margin the point is brought this amount _inside_ the active area
     *        _(default: 0)_
     * @return a projection displacement
     * @see `DeltaFromActivePlane(WidthDepthProjection_t const&, double, double)`
     * 
     * This is the implementation with default values for margins of
     * `DeltaFromActivePlane()`.
     * The depth and width margins are the same, and 0 by default.
     */
    WidthDepthProjection_t DeltaFromActivePlane
      (WidthDepthProjection_t const& proj, double margin = 0.0) const
      { return DeltaFromActivePlane(proj, margin, margin); }
    
    /**
     * @brief Returns the projection, moved onto the plane if necessary.
     * @param proj projection to be checked and moved
     * @return the new value of the projection
     * @see isProjectionOnPlane(), Width(), Height()
     * 
     * The projection proj is defined as in the output of
     * `PointWidthDepthProjection()`.
     * The method caps width and depth of the projection so that it stays on
     * the plane. A new capped value is returned.
     * Since the reference point of the frame is defined as the center of the
     * plane, this action is equivalent to force the width component in
     * @f$ \left[ -\frac{w}{2}, \frac{w}{2} \right] @f$ range and the depth
     * component into @f$ \left[ -\frac{d}{2}, \frac{d}{2} \right] @f$, with
     * @f$ w @f$ and @f$ d @f$ the width and depth of the wire plane.
     */
    WidthDepthProjection_t MoveProjectionToPlane
      (WidthDepthProjection_t const& proj) const;
    
    //@{
    /**
     * @brief Returns the point, moved so that its projection is over the plane.
     * @param point point to be checked and moved
     * @return the new value of the point
     * @see isProjectionOnPlane(), MoveProjectionToPlane(), Width(), Height()
     * 
     * If the projection of the point on the plane falls outside it, the
     * returned point is translated so that its projection is now on the border
     * of the plane. The translation happens along the directions of the plane
     * frame, as described in MoveProjectionToPlane().
     */
    geo::Point_t MovePointOverPlane(geo::Point_t const& point) const;
    TVector3 MovePointOverPlane(TVector3 const& point) const;
    //@}
    
    /**
     * @brief Returns the 3D vector from composition of projection and distance.
     * @tparam Point type of point to be produced (current default is `TVector3`)
     * @param decomp decomposed point
     * @return the 3D vector from composition of projection and distance
     * @see DecomposePointWidthDepth(),
     *      ComposePointWidthDepth(double, DecomposedVector_t::Projection_t const&)
     * 
     * See
     * `ComposePointWidthDepth(double, DecomposedVector_t::Projection_t const&)`
     * for details.
     */
    template <typename Point = DefaultPoint_t>
    Point ComposePoint(WDDecomposedVector_t const& decomp) const
      { return geo::vect::convertTo<Point>(fDecompFrame.ComposePoint(decomp)); }
    
    /**
     * @brief Returns the 3D point from composition of projection and distance.
     * @tparam Point type of point to be produced (current default is `TVector3`)
     * @param distance distance of the target point from the wire plane
     * @param proj projection of the target point on the wire plane
     * @return the 3D vector from composition of projection and distance
     * @see DecomposePointWidthDepth()
     * 
     * The returned vector is the sum of two 3D vectors:
     * 
     * 1. a vector parallel to the plane normal, with norm the input distance
     * 2. a vector lying on the plane, whose projection via
     *    `PointWidthDepthProjection()` gives the input projection
     * 
     * Given the arbitrary definition of the projection reference, it is assumed
     * that the same convention is used as in PointWidthDepthProjection() and
     * DecomposePointWidthDepth().
     * 
     */
    template <typename Point = DefaultPoint_t>
    Point ComposePoint
      (double distance, WidthDepthProjection_t const& proj) const
      { return geo::vect::convertTo<Point>(fDecompFrame.ComposePoint(distance, proj)); }
    
    
    /// @}
    
    
    /// @{
    /**
     * @name Coordinate transformation
     * 
     * Local points and displacement vectors are described by the types
     * `geo::PlaneGeo::LocalPoint_t` and `geo::PlaneGeo::LocalVector_t`,
     * respectively.
     */
    
    /// Transform point from local plane frame to world frame.
    void LocalToWorld(const double* plane, double* world) const
      { fTrans.LocalToWorld(plane, world); }
      
    /// Transform point from local plane frame to world frame.
    TVector3 LocalToWorld(const TVector3& local) const
      { return fTrans.LocalToWorld<TVector3>(local); }
    
    /// Transform point from local plane frame to world frame.
    geo::Point_t toWorldCoords(LocalPoint_t const& local) const
      { return fTrans.toWorldCoords(local); }
    
    /// Transform direction vector from local to world.
    void LocalToWorldVect(const double* plane, double* world) const
      { fTrans.LocalToWorldVect(plane, world); }
    
    /// Transform direction vector from local to world.
    TVector3 LocalToWorldVect(const TVector3& local) const
      { return fTrans.LocalToWorldVect<TVector3>(local); }
    
    /// Transform direction vector from local to world.
    geo::Vector_t toWorldCoords(LocalVector_t const& local) const
      { return fTrans.toWorldCoords(local); }
    
    /// Transform point from world frame to local plane frame.
    void WorldToLocal(const double* world, double* plane) const
      { fTrans.WorldToLocal(world, plane); }
    
    /// Transform point from world frame to local plane frame.
    TVector3 WorldToLocal(TVector3 const& world) const
      { return fTrans.WorldToLocal<TVector3>(world); }
    
    /// Transform point from world frame to local plane frame.
    LocalPoint_t toLocalCoords(geo::Point_t const& world) const
      { return fTrans.toLocalCoords(world); }
    
    /// Transform direction vector from world to local.
    void WorldToLocalVect(const double* world, double* plane) const
      { fTrans.WorldToLocalVect(world, plane); }
    
    /// Transform direction vector from world to local.
    TVector3 WorldToLocalVect(TVector3 const& world) const
      { return fTrans.WorldToLocalVect<TVector3>(world); }
    
    /// Transform direction vector from world to local.
    LocalVector_t toLocalCoords(geo::Vector_t const& world) const
      { return fTrans.toLocalCoords(world); }
    
    /// @}
    
    
    
    /// @{
    /// @name Setters
    
    /// Set the signal view (for TPCGeo).
    void SetView(geo::View_t view)                                  { fView = view; }
    
    /// @}
    
    /// Apply sorting to WireGeo objects.
    void SortWires(geo::GeoObjectSorter const& sorter);
    
    /// Performs all needed updates after the TPC has sorted the planes.
    void UpdateAfterSorting
      (geo::PlaneID planeid, geo::BoxBoundedGeo const& TPCbox);
    
    /// Returns the name of the specified view.
    static std::string ViewName(geo::View_t view);
    
    /// Returns the name of the specified orientation.
    static std::string OrientationName(geo::Orient_t orientation);
    
    
  private:
    
    void FindWire(GeoNodePath_t& path, size_t depth);
    void MakeWire(GeoNodePath_t& path, size_t depth);
    
    /// Sets the geometry directions.
    void DetectGeometryDirections();
    
    /// Returns a direction normal to the plane (pointing is not defined).
    geo::Vector_t GetNormalAxis() const;
    
    /// Updates the cached normal to plane versor; needs the TPC box coordinates.
    void UpdatePlaneNormal(geo::BoxBoundedGeo const& TPCbox);
    
    /// Updates the cached depth and width direction.
    void UpdateWidthDepthDir();
    
    /// Updates the cached direction to increasing wires.
    void UpdateIncreasingWireDir();
    
    /// Updates the cached direction to wire.
    void UpdateWireDir();
    
    /// Updates plane orientation.
    void UpdateOrientation();
    
    /// Updates the stored wire pitch.
    void UpdateWirePitch();
    
    /// Updates the stored wire plane center.
    void UpdateWirePlaneCenter();
    
    /// Updates the stored @f$ \phi_{z} @f$.
    void UpdatePhiZ();
      
    /// Updates the stored view
    void UpdateView();
    
    /// Updates the stored wire pitch with a slower, more robust algorithm.
    void UpdateWirePitchSlow();
    
    /// Updates the position of the wire coordinate decomposition.
    void UpdateDecompWireOrigin();
    
    /// Updates the internally used active area.
    void UpdateActiveArea();
    
    /// Whether the specified wire should have start and end swapped.
    bool shouldFlipWire(geo::WireGeo const& wire) const;
    
  private:
    using WireCollection_t = std::vector<geo::WireGeo>;
    
    using LocalTransformation_t
      = geo::LocalTransformationGeo<TGeoHMatrix, LocalPoint_t, LocalVector_t>;
    
    struct RectSpecs {
      double halfWidth;
      double halfDepth;
      
      double HalfWidth() const { return halfWidth; }
      double HalfDepth() const { return halfDepth; }
      double Width() const { return 2.0 * HalfWidth(); }
      double Depth() const { return 2.0 * HalfDepth(); }
    }; // RectSpecs
    
    LocalTransformation_t fTrans;       ///< Plane to world transform.
    TGeoVolume const*     fVolume;      ///< Plane volume description.
    View_t                fView;        ///< Does this plane measure U, V, or W?
    Orient_t              fOrientation; ///< Is the plane vertical or horizontal?
    WireCollection_t      fWire;        ///< List of wires in this plane.
    double                fWirePitch;   ///< Pitch of wires in this plane.
    double                fSinPhiZ;     ///< Sine of @f$ \phi_{z} @f$.
    double                fCosPhiZ;     ///< Cosine of @f$ \phi_{z} @f$.
    
    geo::Vector_t         fNormal;      ///< Normal to the plane, inward in TPC.
    /// Decomposition on wire coordinates; the main direction is along the wire,
    /// the secondary one is the one measured by the wire, the normal matches
    /// the plane's normal.
    WireDecomposer_t      fDecompWire;
    /// Decomposition on frame coordinates; the main direction is a "width",
    /// the secondary one is just orthogonal to it ("depth").
    /// Normal can differ in sign from the plane one.
    WidthDepthDecomposer_t fDecompFrame;
    RectSpecs             fFrameSize;   ///< Size of the frame of the plane.
    /// Area covered by wires in frame base.
    Rect                  fActiveArea;
    /// Center of the plane, lying on the wire plane.
    geo::Point_t          fCenter;

    geo::PlaneID          fID;          ///< ID of this plane.
    
    friend details::ActiveAreaCalculator;
    
  }; // class PlaneGeo
  
} // namespace geo


//------------------------------------------------------------------------------
//--- template implementation
//---
template <typename Stream>
void geo::PlaneGeo::PrintPlaneInfo(
  Stream&& out,
  std::string indent /* = "" */,
  unsigned int verbosity /* = 1 */
) const {
  
  //----------------------------------------------------------------------------
  out << "plane " << std::string(ID());
  
  if (verbosity-- <= 0) return; // 0
  
  //----------------------------------------------------------------------------
  out
    << " at " << GetCenter<geo::Vector_t>() << " cm"
    << ", theta: " << ThetaZ() << " rad";
  
  if (verbosity-- <= 0) return; // 1
  
  //----------------------------------------------------------------------------
  unsigned int const nWires = Nwires();
  
  out << "\n" << indent
    << "normal to wire: " << PhiZ() << " rad"
      << ", with orientation " << OrientationName(Orientation())
      << ", has " << nWires << " wires measuring " << ViewName(View())
      << " with a wire pitch of " << WirePitch() << " cm"
    ;
  
  if (verbosity-- <= 0) return; // 2
  
  //----------------------------------------------------------------------------
  auto const& normal = GetNormalDirection<geo::Vector_t>();
  auto const& incrZdir = GetIncreasingWireDirection<geo::Vector_t>();
  auto const& wireNormalDir = fDecompWire.NormalDir();
  out << "\n" << indent
    << "normal to plane: " << normal
    << ", direction of increasing wire number: " << incrZdir
    << " [wire frame normal: " << wireNormalDir << "]"
    << " (" << (WireIDincreasesWithZ()? "increases": "decreases") << " with z)";
    
  if (verbosity-- <= 0) return; // 3
  
  //----------------------------------------------------------------------------
  
  auto const& wireDir = GetWireDirection<geo::Vector_t>();
  auto const& widthDir = WidthDir<geo::Vector_t>();
  auto const& depthDir = DepthDir<geo::Vector_t>();
  auto const& frameNormalDir = fDecompFrame.NormalDir();
  
  out << "\n" << indent
    << "wire direction: " << wireDir
    << "; width " << Width() << " cm in direction: " << widthDir
    << ", depth " << Depth() << " cm in direction: " << depthDir
    << " [normal: " << frameNormalDir << "]"
    ;
    
  if (verbosity-- <= 0) return; // 4
  
  //----------------------------------------------------------------------------
  // get the area spanned by the wires
  out << "\n" << indent << "wires cover width "
    << ActiveArea().width.lower << " to " << ActiveArea().width.upper
    << ", depth "
    << ActiveArea().depth.lower << " to " << ActiveArea().depth.upper
    << " cm";
  if (verbosity-- <= 0) return; // 5
  
  //----------------------------------------------------------------------------
  // print also the containing box
  auto const box = BoundingBox();
  out << "\n" << indent
    << "bounding box: " << box.Min() << " -- " << box.Max();
  
//  if (verbosity-- <= 0) return; // 6
  
  //----------------------------------------------------------------------------
} // geo::PlaneGeo::PrintPlaneInfo()


#endif // LARCOREALG_GEOMETRY_PLANEGEO_H
////////////////////////////////////////////////////////////////////////
