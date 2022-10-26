/**
 * @file   larcorealg/Geometry/GeometryCore.h
 * @brief  Access the description of detector geometry
 * @author brebel@fnal.gov
 * @see    larcorealg/Geometry/GeometryCore.cxx
 * @ingroup Geometry
 *
 * Structure of the header:
 *
 *   namespace geo {
 *
 *     // GeometryData_t definition (part of GeometryCore)
 *
 *     // GeometryCore declaration
 *   }
 *
 */

#ifndef LARCOREALG_GEOMETRY_GEOMETRYCORE_H
#define LARCOREALG_GEOMETRY_GEOMETRYCORE_H

// LArSoft libraries
#include "larcorealg/CoreUtils/RealComparisons.h"
#include "larcorealg/CoreUtils/span.h"
#include "larcorealg/Geometry/AuxDetGeo.h"
#include "larcorealg/Geometry/AuxDetSensitiveGeo.h"
#include "larcorealg/Geometry/BoxBoundedGeo.h"
#include "larcorealg/Geometry/ChannelMapAlg.h"
#include "larcorealg/Geometry/CryostatGeo.h"
#include "larcorealg/Geometry/GeoObjectSorter.h"
#include "larcorealg/Geometry/GeometryBuilder.h"
#include "larcorealg/Geometry/GeometryData.h"
#include "larcorealg/Geometry/GeometryDataContainers.h" // geo::TPCDataContainer
#include "larcorealg/Geometry/OpDetGeo.h"
#include "larcorealg/Geometry/PlaneGeo.h"
#include "larcorealg/Geometry/ReadoutDataContainers.h" // readout::ROPDataContainer
#include "larcorealg/Geometry/TPCGeo.h"
#include "larcorealg/Geometry/WireGeo.h"
#include "larcorealg/Geometry/details/geometry_iterators.h"
#include "larcorealg/Geometry/fwd.h"
#include "larcorealg/Geometry/geo_vectors_utils.h"       // geo::vect namespace
#include "larcoreobj/SimpleTypesAndConstants/RawTypes.h" // raw::ChannelID_t
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"
#include "larcoreobj/SimpleTypesAndConstants/readout_types.h"

// Framework and infrastructure libraries
#include "fhiclcpp/ParameterSet.h"

// ROOT libraries
#include "TVector3.h"

// C/C++ standard libraries
#include <cstddef>  // size_t
#include <iterator> // std::forward_iterator_tag
#include <memory>   // std::shared_ptr<>
#include <set>
#include <string>
#include <type_traits> // std::is_base_of<>
#include <utility>
#include <vector>

// ROOT class prototypes
class TGeoManager;
class TGeoNode;
class TGeoVolume;
class TGeoMaterial;

namespace unrelated {
  // Functions to allow determination if two wires intersect, and if so where.
  // This is useful information during 3D reconstruction.
  //......................................................................
  bool ValueInRange(double value, double min, double max);

  // The following functions are utilized to determine if two wires
  // in the TPC intersect or not, and if they do then
  // determine the coordinates of the intersection.

  /**
     * @brief Computes the intersection between two lines on a plane
     * @param A_start_x x coordinate of one point of the first segment
     * @param A_start_y y coordinate of one point of the first segment
     * @param A_end_x x coordinate of another point of the first segment
     * @param A_end_y y coordinate of another point of the first segment
     * @param B_start_x x coordinate of one point of the second segment
     * @param B_start_y y coordinate of one point of the second segment
     * @param B_end_x x coordinate of another point of the second segment
     * @param B_end_y y coordinate of another point of the second segment
     * @param x _(output)_ variable to store the x coordinate of intersection
     * @param y _(output)_ variable to store the y coordinate of intersection
     * @return whether intersection exists
     *
     * The order of the ends is not relevant.
     * The return value is `false` if the two segments are parallel.
     * In that case, `x` and `y` variables are not changed.
     * Otherwise, they hold the intersection coordinate, even if the
     * intersection point is beyond one or both the segments.
     */
  bool IntersectLines(double A_start_x,
                      double A_start_y,
                      double A_end_x,
                      double A_end_y,
                      double B_start_x,
                      double B_start_y,
                      double B_end_x,
                      double B_end_y,
                      double& x,
                      double& y);

  /**
     * @brief Computes the intersection between two segments on a plane
     * @param A_start_x x coordinate of the start of the first segment
     * @param A_start_y y coordinate of the start of the first segment
     * @param A_end_x x coordinate of the end of the first segment
     * @param A_end_y y coordinate of the end of the first segment
     * @param B_start_x x coordinate of the start of the second segment
     * @param B_start_y y coordinate of the start of the second segment
     * @param B_end_x x coordinate of the end of the second segment
     * @param B_end_y y coordinate of the end of the second segment
     * @param x _(output)_ variable to store the x coordinate of intersection
     * @param y _(output)_ variable to store the y coordinate of intersection
     * @return whether intersection exists and is on both segments
     *
     * The order of the ends is not relevant.
     * The return value is `false` if the two segments are parallel, or if their
     * intersection point is not on _both_ the segments.
     * If the segments are parallel, x and y variables are not changed.
     * Otherwise, they hold the intersection coordinate, even if the
     * intersection point is beyond one or both the segments.
     */
  bool IntersectSegments(double A_start_x,
                         double A_start_y,
                         double A_end_x,
                         double A_end_y,
                         double B_start_x,
                         double B_start_y,
                         double B_end_x,
                         double B_end_y,
                         double& x,
                         double& y);
}

/// Namespace collecting geometry-related classes utilities
namespace geo {

  // BEGIN Geometry group ------------------------------------------------------
  /// @ingroup Geometry
  /// @{

  template <typename Iter>
  struct IteratorBox : util::span<Iter> {
    using util::span<Iter>::span;
  };

  template <typename Iter>
  IteratorBox(Iter, Iter)->IteratorBox<Iter>;

  //
  // GeometryCore
  //

  /** **************************************************************************
   * @brief Description of geometry of one entire detector
   *
   * @note All lengths are specified in centimetres
   *
   *
   * How to correctly instantiate a GeometryCore object
   * ---------------------------------------------------
   *
   * Instantiation is a multi-step procedure:
   * 1. construct a GeometryCore object (the "service provider"),
   *    with the full configuration; at this step, configuration is just stored
   * 2. load a geometry with GeometryCore::LoadGeometryFile();
   *    this loads the detector geometry information
   * 3. prepare a channel map algorithm object (might use for example
   *    GeometryCore::DetectorName() or the detector geometry from the
   *    newly created object, but any use of channel mapping related functions
   *    is forbidden and it would yield undefined behaviour (expected to be
   *    catastrophic)
   * 4. acquire the channel mapping algorithm with
   *    GeometryCore::ApplyChannelMap().
   *
   * Step 3 (creation of the channel mapping algorithm object) can be performed
   * at any time before step 4, provided that no GeometryCore instance is needed
   * for it.
   *
   *
   * Configuration parameters
   * -------------------------
   *
   * - *Name* (string; mandatory): string identifying the detector; it can be
   *   different from the base name of the file used to initialize the geometry;
   *   standard names are recommended by each experiment.
   *   This name can be used, for example, to select which channel mapping
   *   algorithm to use.
   * - *SurfaceY* (real; mandatory): depth of the detector, in centimetrs;
   *   see SurfaceY() for details
   * - *MinWireZDist* (real; default: 3)
   * - *PositionEpsilon* (real; default: 0.01%) set the default tolerance
   *   (see DefaultWiggle())
   *
   */
  class GeometryCore {
    using DefaultVector_t = TVector3; ///< Default template argument.
    using DefaultPoint_t = TVector3;  ///< Default template argument.

  public:
    /// Simple class with two points (a pair with aliases).
    template <typename Point>
    struct Segment : public std::pair<Point, Point> {

      // use the base class constructors
      using std::pair<Point, Point>::pair;

      Point const& start() const { return this->first; }
      Point& start() { return this->first; }

      Point const& end() const { return this->second; }
      Point& end() { return this->second; }

    }; // struct Segment_t

    using Segment_t = Segment<DefaultPoint_t>;

    /// Type of list of cryostats
    using CryostatList_t = GeometryData_t::CryostatList_t;
    /// Type of list of auxiliary detectors
    using AuxDetList_t = GeometryData_t::AuxDetList_t;

    /// Wires must be found in GDML description within this number of nested
    /// volumes.
    static constexpr std::size_t MaxWireDepthInGDML = 20U;

    /**
     * @brief Initialize geometry from a given configuration
     * @param pset configuration parameters
     *
     * This constructor does not load any geometry description.
     * The next step is to do exactly that, by GeometryCore::LoadGeometryFile().
     */
    GeometryCore(fhicl::ParameterSet const& pset);

    /// Destructor
    ~GeometryCore();

    // this object is not copiable nor moveable (see also issue #14384);
    // currently, auxiliary detectors are stored as bare pointers,
    // which prevents trivial copy or move.
    GeometryCore(GeometryCore const&) = delete;
    GeometryCore(GeometryCore&&) = delete;
    GeometryCore& operator=(GeometryCore const&) = delete;
    GeometryCore& operator=(GeometryCore&&) = delete;

    /**
     * @brief Returns the tolerance used in looking for positions
     * @return the tolerance value
     *
     * This parameter is used as tolerance ("wiggle") for methods that require
     * it (e.g. `geo::CryostatGeo::FindTPCAtPosition()`).
     * Typically, it's a additional fraction of tolerance: 0 means no tolerance,
     * 0.1 means 10% tolerance.
     *
     * @todo Confirm the definition of wiggle: this one is taken from other doc
     */
    double DefaultWiggle() const { return fPositionWiggle; }

    /**
     * @brief Returns the full directory path to the geometry file source
     * @return the full directory path to the geometry file source
     *
     * This is the full path of the source of the detector geometry GeometryCore
     * relies on.
     */
    std::string const& ROOTFile() const { return fROOTfile; }

    /**
     * @brief Returns the full directory path to the GDML file source
     * @return the full directory path to the GDML file source
     *
     * This is the full path of the source of the detector geometry handed to
     * the detector simulation (GEANT).
     */
    std::string const& GDMLFile() const { return fGDMLfile; }

    // BEGIN Detector information
    /// @name Detector information
    /// @{

    //
    // global features
    //
    /// Returns a string with the name of the detector, as configured
    std::string const& DetectorName() const { return fDetectorName; }

    //
    // position
    //

    /// Returns a pointer to the world volume.
    TGeoVolume const* WorldVolume() const;

    /**
     * @brief Fills the arguments with the boundaries of the world
     * @param xlo (output) pointer to the lower x coordinate
     * @param xlo (output) pointer to the upper x coordinate
     * @param ylo (output) pointer to the lower y coordinate
     * @param ylo (output) pointer to the upper y coordinate
     * @param zlo (output) pointer to the lower z coordinate
     * @param zlo (output) pointer to the upper z coordinate
     * @throw cet::exception (`"GeometryCore"` category) if no world found
     * @see `GetWorldVolumeName()`
     *
     * This method fills the boundaries of the world volume
     * (`GetWorldVolumeName()`).
     *
     * If a pointer is null, its coordinate is skipped.
     *
     * @deprecated Use the version without arguments instead.
     */
    void WorldBox(double* xlo, double* xhi, double* ylo, double* yhi, double* zlo, double* zhi)
      const;

    /// Returns a box with the extremes of the world volume (from shape axes).
    /// @see `GetWorldVolumeName()`
    BoxBoundedGeo WorldBox() const;

    /**
     * @brief The position of the detector respect to earth surface
     * @return typical y position at surface in units of cm
     *
     * This is the depth (y) of the surface (where earth meets air) for this
     * detector site.
     * The number is expressed in world coordinates and in centimetres,
     * and it represents the y coordinate of earth surface.
     * A negative value means that the origin of coordinates, typically matching
     * the detector centre, is above surface.
     *
     * @todo check that this is actually how it is used
     */
    //
    Length_t SurfaceY() const { return fSurfaceY; }

    //
    // object description and information
    //

    /// Access to the ROOT geometry description manager
    TGeoManager* ROOTGeoManager() const;

    /// Return the name of the world volume (needed by Geant4 simulation)
    const std::string GetWorldVolumeName() const;

    /// Returns the absolute  coordinates of the detector enclosure volume [cm].
    /// @param name name of the volume to be sought (default: `volDetEnclosure`)
    /// @throw cet::exception if the specified volume is not found
    BoxBoundedGeo DetectorEnclosureBox(std::string const& name = "volDetEnclosure") const;

    //@{
    /**
     * @brief Returns the name of the deepest volume containing specified point
     * @param point the location to query, in world coordinates
     * @return name of the volume containing the point
     *
     * @todo what happens if none?
     * @todo Unify the coordinates type
     */
    std::string VolumeName(Point_t const& point) const;
    std::string VolumeName(TVector3 const& point) const { return VolumeName(vect::toPoint(point)); }
    //@}

    /**
     * @brief Returns all the nodes with volumes with any of the specified names
     * @param vol_names list of names of volumes
     * @return list of nodes found
     *
     * All the nodes in the geometry are checked, and all the ones that contain
     * a volume with a name among the ones specified in vol_names are saved
     * in the collection and returned.
     */
    std::vector<TGeoNode const*> FindAllVolumes(std::set<std::string> const& vol_names) const;

    /**
     * @brief Returns paths of all nodes with volumes with the specified names
     * @param vol_names list of names of volumes
     * @return list paths of the found nodes
     *
     * All the nodes in the geometry are checked, and the path of all the ones
     * that contain a volume with a name among the ones specified in vol_names
     * is saved in the collection and returned.
     * A node path is a ordered list of all nodes leading to the final one,
     * starting from thetop level (root) down. The node at the `back()` of the
     * path is the one with name in vol_names.
     * No empty paths are returned.
     */
    std::vector<std::vector<TGeoNode const*>> FindAllVolumePaths(
      std::set<std::string> const& vol_names) const;

    /// Returns the material at the specified position
    TGeoMaterial const* Material(Point_t const& point) const;
    //@{
    /**
     * @brief Name of the deepest material containing the point xyz
     * @return material of the origin by default
     */
    std::string MaterialName(TVector3 const& point) const
    {
      return MaterialName(vect::toPoint(point));
    }
    std::string MaterialName(Point_t const& point) const;
    //@}

    //@{
    /// Returns the total mass [kg] of the specified volume (default: world).
    double TotalMass() const { return TotalMass(GetWorldVolumeName()); }
    double TotalMass(std::string vol) const;
    //@}

    //@{
    /**
     * @brief Returns the column density between two points.
     * @param p1 the first point
     * @param p2 the second point
     * @return the column density [kg / cm&sup2;]
     *
     * The column density is defined as
     * @f$ \int_{\vec{p}_{1}}^{\vec{p}_{2}} \rho(\vec{p}) d\vec{p} @f$
     * where @f$ \rho(\vec{p}) @f$ is the density at point @f$ \vec{p} @f$,
     * which the integral leads from `p1` to `p2` in a straight line.
     *
     * Both points are specified in world coordinates.
     */
    double MassBetweenPoints(Point_t const& p1, Point_t const& p2) const;
    double MassBetweenPoints(double* p1, double* p2) const;
    //@}

    /// Prints geometry information with maximum verbosity.
    template <typename Stream>
    void Print(Stream&& out, std::string indent = "  ") const;

    /// @brief Returns a string with complete geometry information.
    /// @see `Print()`
    std::string Info(std::string indent = "  ") const;

    /// @}
    // END Detector information

    /**
     * @brief Returns the ID of the first element of the detector.
     * @tparam GeoID type of the ID to be returned
     * @return ID of the first subelement in the detector
     */
    template <typename GeoID>
    GeoID GetBeginID() const
    {
      GeoID id;
      GetBeginID(id);
      return id;
    }

    /**
     * @brief Returns the ID next to the specified one.
     * @tparam GeoID type of the ID to be returned
     * @param id the element ID to be incremented
     * @return ID of the next subelement after `id`
     */
    template <typename GeoID>
    GeoID GetNextID(GeoID const& id) const
    {
      auto nextID(id);
      IncrementID(nextID);
      return nextID;
    }

    /**
     * @brief Returns the (possibly invalid) ID after the last subelement of
     *        the detector.
     * @tparam GeoID type of the ID to be returned
     * @return ID after the last subelement in the specified geometry element
     */
    template <typename GeoID>
    GeoID GetEndID() const
    {
      GeoID id;
      GetEndID(id);
      return id;
    }

    /**
     * @brief Returns the ID of the first subelement of the specified element.
     * @tparam GeoID type of the ID to be returned
     * @tparam ContextID type of the ID of the containing element
     * @param id ID of the containing element
     * @return ID of the first subelement in the specified geometry element
     */
    template <typename GeoID, typename ContextID>
    GeoID GetBeginID(ContextID const& id) const;

    /**
     * @brief Returns the (possibly invalid) ID after the last subelement of
     *        the specified element.
     * @tparam GeoID type of the ID to be returned
     * @tparam ContextID type of the ID of the containing element
     * @param id ID of the containing element
     * @return ID (possibly invalid) after the last subelement in the
     *         specified geometry element
     */
    template <typename GeoID, typename ContextID>
    GeoID GetEndID(ContextID const& id) const;

    /// @name Cryostat access and information
    /// @{

    //
    // group features
    //

    //@{
    /**
     * @brief Returns the number of cryostats in the detector
     *
     * The NElements() and NSiblingElements() methods are overloaded and their
     * return depends on the type of ID.
     *
     * @todo Change return type to size_t
     */
    unsigned int Ncryostats() const { return Cryostats().size(); }
    unsigned int NElements() const { return Ncryostats(); }
    unsigned int NSiblingElements(CryostatID const&) const { return Ncryostats(); }
    //@}

    //
    // access
    //

    //@{
    /**
     * @brief Returns whether we have the specified cryostat
     *
     * The HasElement() method is overloaded and its meaning depends on the type
     * of ID.
     */
    bool HasCryostat(CryostatID const& cryoid) const { return cryoid.Cryostat < Ncryostats(); }
    bool HasElement(CryostatID const& cryoid) const { return HasCryostat(cryoid); }
    //@}

    //@{
    /**
     * @brief Returns the specified cryostat
     * @param cstat number of cryostat
     * @param cryoid cryostat ID
     * @return a constant reference to the specified cryostat
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     *
     * The GetElement() method is overloaded and its return depends on the type
     * of ID.
     *
     * @todo Make the cryostat number mandatory (as CryostatID)
     */
    CryostatGeo const& Cryostat(CryostatID const& cryoid) const;
    CryostatGeo const& Cryostat(unsigned int const cstat = 0) const
    {
      return Cryostat(CryostatID(cstat));
    }
    CryostatGeo const& GetElement(CryostatID const& cryoid) const { return Cryostat(cryoid); }
    //@}

    //@{
    /**
     * @brief Returns the specified cryostat
     * @param cryoid cryostat ID
     * @return a constant pointer to the specified cryostat, or nullptr if none
     *
     * The GetElementPtr() method is overloaded and its return depends on the
     * type of ID.
     */
    CryostatGeo const* CryostatPtr(CryostatID const& cryoid) const
    {
      return HasCryostat(cryoid) ? &(Cryostats()[cryoid.Cryostat]) : nullptr;
    }
    CryostatGeo const* GetElementPtr(CryostatID const& cryoid) const { return CryostatPtr(cryoid); }
    //@}

    /**
     * @brief Returns the cryostat at specified location.
     * @param point the location [cm]
     * @return pointer to the `geo::CryostatGeo` including `point`, or `nullptr`
     *
     * The tolerance used here is the one returned by DefaultWiggle().
     */
    CryostatGeo const* PositionToCryostatPtr(Point_t const& point) const;

    /**
     * @brief Returns the ID of the cryostat at specified location.
     * @param point the location [cm]
     * @return ID of the cryostat including `point` (invalid if none)
     *
     * The tolerance used here is the one returned by DefaultWiggle().
     */
    CryostatID PositionToCryostatID(Point_t const& point) const;

    //@{
    /**
     * @brief Returns the cryostat at specified location.
     * @param point the location [cm]
     * @return a constant reference to the `geo::CryostatGeo` containing `point`
     * @throws cet::exception ("Geometry" category) if no cryostat matches
     *
     * The tolerance used here is the one returned by DefaultWiggle().
     */
    CryostatGeo const& PositionToCryostat(Point_t const& point) const;

    //
    // iterators
    //

    /// Initializes the specified ID with the ID of the first cryostat
    void GetBeginID(CryostatID& id) const { id = CryostatID(0, HasCryostat(CryostatID(0))); }

    /// Initializes the specified ID with the invalid ID after the last cryostat
    void GetEndID(CryostatID& id) const { id = CryostatID(Ncryostats(), false); }

    /// Sets the ID to the ID after the specified one.
    /// @return whether the ID is actually valid (validity flag is also set)
    bool IncrementID(CryostatID& id) const; // inline implementation

    /// Returns an iterator pointing to the first cryostat ID

    template <typename BaseID, typename GeoID>
    static constexpr bool is_base_of_strict{std::is_base_of<BaseID, GeoID>{} &&
                                            !std::is_same<BaseID, GeoID>{}};

    template <typename T>
    auto begin() const
    {
      using namespace details;
      if constexpr (std::is_base_of<CryostatID, T>{}) {
        return id_iterator<T>{this, GetBeginID<T>()};
      }
      else {
        using ID_t = typename T::ID_t;
        return element_iterator_for<T>{this, id_iterator<ID_t>{this, GetBeginID<ID_t>()}};
      }
    }

    template <typename T>
    auto end() const
    {
      using namespace details;
      if constexpr (std::is_base_of<CryostatID, T>{}) {
        return id_iterator<T>{this, GetEndID<T>()};
      }
      else {
        using ID_t = typename T::ID_t;
        return element_iterator_for<T>{this, id_iterator<ID_t>{this, GetEndID<ID_t>()}};
      }
    }

    template <typename T, typename BaseID>
    auto begin(BaseID const& id) const
    {
      using namespace details;
      if constexpr (std::is_base_of<CryostatID, T>{}) {
        static_assert(is_base_of_strict<BaseID, T>);
        return id_iterator<T>{this, GetBeginID<T>(id)};
      }
      else {
        using ID_t = typename T::ID_t;
        static_assert(is_base_of_strict<BaseID, ID_t>);
        return element_iterator_for<T>{this, id_iterator<ID_t>{this, GetBeginID<ID_t>(id)}};
      }
    }

    template <typename T, typename BaseID>
    auto end(BaseID const& id) const
    {
      using namespace details;
      if constexpr (std::is_base_of<CryostatID, T>{}) {
        static_assert(is_base_of_strict<BaseID, T>);
        return id_iterator<T>{this, GetEndID<T>(id)};
      }
      else {
        using ID_t = typename T::ID_t;
        static_assert(is_base_of_strict<BaseID, ID_t>);
        return element_iterator_for<T>{this, id_iterator<ID_t>{this, GetEndID<ID_t>(id)}};
      }
    }

    // template <typename ID, typename BaseID>
    // details::id_iterator<ID> begin_id(BaseID const& base_id) const
    // {
    //   return {this, GetBeginID<ID>(base_id)};
    // }

    // template <typename ID, typename BaseID>
    // details::id_iterator<ID> end_id(BaseID const& base_id) const
    // {
    //   static_assert(is_base_of_strict<BaseID, ID>);
    //   return {this, GetEndID<ID>(base_id)};
    // }

    // template <typename Element, typename BaseID>
    // details::element_iterator_for<Element> begin(BaseID const& base_id) const
    // {
    //   return {this, begin_id<typename Element::ID_t>(base_id)};
    // }

    // template <typename Element, typename BaseID>
    // details::element_iterator_for<Element> end(BaseID const& base_id) const
    // {
    //   return {this, end_id<typename Element::ID_t>(base_id)};
    // }

    template <typename T>
    auto Iterate() const
    {
      return IteratorBox{begin<T>(), end<T>()};
    }

    template <typename T, typename ID>
    auto Iterate(ID const& id) const
    {
      return IteratorBox{begin<T>(id), end<T>(id)};
    }

    //
    // single object features
    //

    //@{
    /// Returns the half width of the cryostat (x direction)
    Length_t CryostatHalfWidth(CryostatID const& cid) const;
    Length_t CryostatHalfWidth(unsigned int cstat = 0) const
    {
      return CryostatHalfWidth(CryostatID(cstat));
    }
    //@}

    //@{
    /// Returns the height of the cryostat (y direction)
    Length_t CryostatHalfHeight(CryostatID const& cid) const;
    Length_t CryostatHalfHeight(unsigned int cstat = 0) const
    {
      return CryostatHalfHeight(CryostatID(cstat));
    }
    //@}

    //@{
    /// Returns the length of the cryostat (z direction)
    Length_t CryostatLength(CryostatID const& cid) const;
    Length_t CryostatLength(unsigned int cstat = 0) const
    {
      return CryostatLength(CryostatID(cstat));
    }
    //@}

    //
    // object description
    //

    //@{
    /**
     * @brief Return the name of LAr TPC volume
     * @param cstat index of the cryostat
     * @return the name of the specified TPC
     *
     * This information is used in the event display.
     *
     * @todo Use a cryostat ID instead
     * @todo What if it does not exist?
     */
    std::string GetCryostatVolumeName(CryostatID const& cid) const;
    std::string GetCryostatVolumeName(unsigned int const cstat = 0) const
    {
      return GetCryostatVolumeName(CryostatID(cstat));
    }
    //@}

    /// @} Cryostat access and information

    /// @name TPC access and information
    /// @{

    //
    // group features
    //

    /**
     * @brief Returns the total number of TPCs in the specified cryostat
     * @param cstat cryostat number
     *
     * @todo Make the cryostat number mandatory (as CryostatID)
     * @todo Change return type to size_t
     * @todo what happens if it does not exist?
     */
    unsigned int NTPC(unsigned int cstat = 0) const { return NTPC(CryostatID(cstat)); }

    /// Returns the largest number of TPCs a cryostat in the detector has
    unsigned int MaxTPCs() const;

    /// Returns the total number of TPCs in the detector
    unsigned int TotalNTPC() const;

    /**
     * @brief Returns a container with one entry per TPC.
     * @tparam T type of data in the container
     * @return a container with one default-constructed `T` per TPC
     * @see `geo::TPCDataContainer`
     *
     * The working assumption is that all cryostats have the same number of
     * TPCs. It is always guaranteed that all existing TPCs have an entry in
     * the container, although if the previous working assumption is not
     * satisfied there will be entries in the containers which are not
     * associated to a valid TPC.
     *
     * The interface of the container is detailed in the documentation of the
     * container itself, `geo::TPCDataContainer`. Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto tracksPerTPC
     *   = geom->makeTPCData<std::vector<recob::Track const*>>();
     *
     * for (recob::Track const& track: tracks) {
     *   geo::TPCGeo const* tpc = geom->PositionToTPCptr(track.Start());
     *   if (tpc) tracksPerTPC[tpc->ID()].push_back(&track);
     * } // for
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * where the container will be filled with pointers to all tracks starting
     * from a given TPC (tracks reconstructed as starting outside the TPCs will
     * be not saved in the container).
     */
    template <typename T>
    TPCDataContainer<T> makeTPCData() const
    {
      return {Ncryostats(), MaxTPCs()};
    }

    /**
     * @brief Returns a container with one entry per TPC.
     * @tparam T type of data in the container
     * @param defValue the initial value of all elements in the container
     * @return a container with a value `defValue` per each TPC
     * @see `geo::TPCDataContainer`
     *
     * This function operates as `makeTPCData() const`, except that copies
     * the specified value into all the entries of the container. Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto nTracksPerTPC = geom->makeTPCData(0U);
     *
     * for (recob::Track const& track: tracks) {
     *   geo::TPCGeo const* tpc = geom->PositionToTPCptr(track.Start());
     *   if (tpc) ++(tracksPerTPC[tpc->ID()]);
     * } // for
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template <typename T>
    TPCDataContainer<T> makeTPCData(T const& defValue) const
    {
      return {Ncryostats(), MaxTPCs(), defValue};
    }

    //@{
    /**
     * @brief Returns the total number of TPCs in the specified cryostat
     * @param cryoid cryostat number
     * @return number of TPCs in specified cryostat, or 0 if no cryostat found
     *
     * The NElements() and NSiblingElements() methods are overloaded and their
     * return depends on the type of ID.
     *
     * @todo Change return type to size_t
     */
    unsigned int NTPC(CryostatID const& cryoid) const
    {
      CryostatGeo const* pCryo = GetElementPtr(cryoid);
      return pCryo ? pCryo->NElements() : 0;
    }
    unsigned int NElements(CryostatID const& cryoid) const { return NTPC(cryoid); }
    unsigned int NSiblingElements(TPCID const& tpcid) const { return NTPC(tpcid); }
    //@}

    //
    // access
    //
    /// Returns whether we have the specified TPC
    bool HasTPC(TPCID const& tpcid) const
    {
      CryostatGeo const* pCryo = CryostatPtr(tpcid);
      return pCryo ? pCryo->HasTPC(tpcid) : false;
    }

    /// Returns whether we have the specified TPC
    bool HasElement(TPCID const& tpcid) const { return HasTPC(tpcid); }

    //@{
    /**
     * @brief Returns the specified TPC
     * @param tpcid ID of the tpc
     * @param tpc tpc number within the cryostat
     * @param cstat number of cryostat
     * @return a constant reference to the specified TPC
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     * @throw cet::exception (`TPCOutOfRange` category) if no such TPC
     *
     * The GetElement() method is overloaded and its return depends on the type
     * of ID.
     *
     * @todo remove the version with integers
     */
    TPCGeo const& TPC(unsigned int const tpc = 0, unsigned int const cstat = 0) const
    {
      return TPC(TPCID(cstat, tpc));
    }
    TPCGeo const& TPC(TPCID const& tpcid) const { return Cryostat(tpcid).TPC(tpcid); }
    TPCGeo const& GetElement(TPCID const& tpcid) const { return TPC(tpcid); }
    //@}

    //@{
    /**
     * @brief Returns the specified TPC
     * @param tpcid TPC ID
     * @return a constant pointer to the specified TPC, or nullptr if none
     *
     * The GetElementPtr() method is overloaded and its return depends on the
     * type of ID.
     */
    TPCGeo const* TPCPtr(TPCID const& tpcid) const
    {
      CryostatGeo const* pCryo = CryostatPtr(tpcid);
      return pCryo ? pCryo->TPCPtr(tpcid) : nullptr;
    }
    TPCGeo const* GetElementPtr(TPCID const& tpcid) const { return TPCPtr(tpcid); }
    //@}

    /**
     * @brief Returns the ID of the TPC at specified location.
     * @param worldLoc 3D coordinates of the point (world reference frame) [cm]
     * @return the TPC ID, or an invalid one if no TPC is there
     */
    TPCID FindTPCAtPosition(double const worldLoc[3]) const
    {
      return FindTPCAtPosition(vect::makePointFromCoords(worldLoc));
    }

    //@{
    /**
     * @brief Returns the ID of the TPC at specified location.
     * @param worldLoc 3D point (world reference frame, centimeters)
     * @return the TPC ID, or an invalid one if no TPC is there
     */
    TPCID FindTPCAtPosition(Point_t const& point) const;
    TPCID FindTPCAtPosition(TVector3 const& point) const
    {
      return FindTPCAtPosition(vect::toPoint(point));
    }
    //@}

    /**
     * @brief Returns the TPC at specified location.
     * @param point the location [cm]
     * @return the `geo::TPCGeo` including `point`, or `nullptr` if none
     */
    TPCGeo const* PositionToTPCptr(Point_t const& point) const;

    //@{
    /**
     * @brief Returns the TPC at specified location.
     * @param point the location [cm]
     * @return a constant reference to the `geo::TPCGeo` including `point`
     * @throws cet::exception ("Geometry" category) if no TPC matches
     */
    TPCGeo const& PositionToTPC(Point_t const& point) const;
    TPCGeo const& PositionToTPC(double const point[3]) const
    {
      return PositionToTPC(vect::makePointFromCoords(point));
    }
    //@}

    /**
     * @brief Returns the ID of the TPC at specified location.
     * @param point the location [cm]
     * @return ID of the TPC at specified location, invalid if none
     * @see `PositionToTPC()`
     */
    TPCID PositionToTPCID(Point_t const& point) const;

    ///
    /// iterators
    ///

    /// Initializes the specified ID with the ID of the first TPC.
    void GetBeginID(TPCID& id) const
    {
      GetBeginID(id.asCryostatID());
      id.TPC = 0;
    }

    /// Initializes the specified ID with the invalid ID after the last TPC.
    void GetEndID(TPCID& id) const;

    /// Sets the ID to the ID after the specified one.
    /// @return whether the ID is actually valid (validity flag is also set)
    bool IncrementID(TPCID& id) const; // inline implementation

    /// Returns the ID of the first TPC in the specified cryostat.
    TPCID GetBeginTPCID(CryostatID const& id) const { return {id, 0}; }

    /// Returns the (possibly invalid) ID after the last TPC of the specified
    /// cryostat.
    TPCID GetEndTPCID(CryostatID const& id) const;

    //
    // single object features
    //

    //@{
    /**
     * @brief Returns the half width of the active volume of the specified TPC.
     * @param tpcid ID of the TPC
     * @param tpc TPC number within the cryostat
     * @param cstat number of cryostat
     * @return the value of the half width of the specified TPC
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     * @throw cet::exception (`TPCOutOfRange` category) if no such TPC
     * @see geo::TPCGeo::ActiveHalfWidth()
     *
     * @todo deprecate this function
     * @todo rename the function
     */
    Length_t DetHalfWidth(TPCID const& tpcid) const;
    Length_t DetHalfWidth(unsigned int tpc = 0, unsigned int cstat = 0) const
    {
      return DetHalfWidth(TPCID(cstat, tpc));
    }
    //@}

    //@{
    /**
     * @brief Returns the half height of the active volume of the specified TPC.
     * @param tpcid ID of the TPC
     * @param tpc TPC number within the cryostat
     * @param cstat number of cryostat
     * @return the value of the half height of the specified TPC
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     * @throw cet::exception (`TPCOutOfRange` category) if no such TPC
     * @see geo::TPCGeo::ActiveHalfHeight()
     *
     * See `geo::TPCGeo::ActiveHalfHeight()` for more details.
     *
     * @todo deprecate this function
     * @todo rename the function
     */
    Length_t DetHalfHeight(TPCID const& tpcid) const;
    Length_t DetHalfHeight(unsigned int tpc = 0, unsigned int cstat = 0) const
    {
      return DetHalfHeight(TPCID(cstat, tpc));
    }
    //@}

    //@{
    /**
     * @brief Returns the length of the active volume of the specified TPC.
     * @param tpcid ID of the TPC
     * @param tpc TPC number within the cryostat
     * @param cstat number of cryostat
     * @return the value of the length of the specified TPC
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     * @throw cet::exception (`TPCOutOfRange` category) if no such TPC
     * @see geo::TPCGeo::ActiveLength()
     *
     * See `geo::TPCGeo::ActiveLength()` for more details.
     *
     * @todo deprecate this function
     * @todo rename the function
     */
    Length_t DetLength(TPCID const& tpcid) const;
    Length_t DetLength(unsigned int tpc = 0, unsigned int cstat = 0) const
    {
      return DetLength(TPCID(cstat, tpc));
    }
    //@}

    //@{
    /**
     * @brief Returns the center of side of the detector facing the beam.
     * @tparam Point _(default: `DefaultPoint_t`)_ return this point type
     * @param tpcid ID of the TPC
     * @return position of center of TPC face toward the beam,
     *         in world coordinates [cm]
     *
     * Effectively, this is the center of the side of TPC active volume
     * which faces the negative _z_ direction, the first that a beam following
     * the positive _z_ direction crosses.
     */
    template <typename Point>
    Point GetTPCFrontFaceCenter(TPCID const& tpcid) const
    {
      return TPC(tpcid).GetFrontFaceCenter<Point>();
    }
    DefaultPoint_t GetTPCFrontFaceCenter(TPCID const& tpcid) const
    {
      return GetTPCFrontFaceCenter<DefaultPoint_t>(tpcid);
    }
    //@}

    //@{
    /**
     * @brief Returns the center of side of the detector facing the beam.
     * @tparam Point _(default: `DefaultPoint_t`)_ return this point type
     * @param tpc _(default: `0`)_ TPC number within the cryostat `cstat`
     * @param cstat _(default: `0`)_ number of cryostat
     * @return position of center of TPC face toward the beam,
     *         in world coordinates [cm]
     * @see `GetTPCFrontFaceCenter(geo::TPCID const&)`
     *
     * @note Please use `GetTPCFrontFaceCenter(geo::TPCID const&)` instead.
     */
    template <typename Point>
    Point GetTPCFrontFaceCenter(unsigned int tpc = 0, unsigned int cstat = 0) const
    {
      return GetTPCFrontFaceCenter<Point>(TPCID(cstat, tpc));
    }
    DefaultPoint_t GetTPCFrontFaceCenter(unsigned int tpc = 0, unsigned int cstat = 0) const
    {
      return GetTPCFrontFaceCenter<DefaultPoint_t>(tpc, cstat);
    }
    //@}

    //
    // object description
    //

    //@{
    /**
     * @brief Return the name of specified LAr TPC volume
     * @param tpcid ID of the TPC
     * @param tpc index of TPC in the cryostat
     * @param cstat index of the cryostat
     * @return the name of the specified TPC
     *
     * This information is used by Geant4 simulation
     *
     * @todo Use a TPCID instead
     * @todo What if it does not exist?
     */
    std::string GetLArTPCVolumeName(TPCID const& tpcid) const;
    std::string GetLArTPCVolumeName(unsigned int const tpc = 0, unsigned int const cstat = 0) const
    {
      return GetLArTPCVolumeName(TPCID(cstat, tpc));
    }
    //@}

    /// @} TPC access and information

    /// @name Plane access and information
    /// @{

    //
    // group features
    //

    /**
     * @brief Returns the total number of wire planes in the specified TPC
     * @param tpc tpc number within the cryostat
     * @param cstat cryostat number
     *
     * @todo Make all the arguments mandatory (as TPCID)
     * @todo Change return type to size_t
     * @todo what happens if TPC does not exist?
     */
    unsigned int Nplanes(unsigned int tpc = 0, unsigned int cstat = 0) const
    {
      return Nplanes(TPCID(cstat, tpc));
    }

    /// Returns the largest number of planes among all TPCs in this detector
    unsigned int MaxPlanes() const;

    /**
     * @brief Returns a container with one entry per wire plane.
     * @tparam T type of data in the container
     * @return a container with one default-constructed `T` per plane
     * @see `geo::PlaneDataContainer`
     *
     * The working assumption is that all cryostats have the same number of
     * TPCs, and all TPCs have the same number of planes. It is always
     * guaranteed that all existing planes have an entry in the container,
     * although if the previous working assumption is not satisfied there will
     * be entries in the containers which are not associated to a valid plane.
     *
     * The interface of the container is detailed in the documentation of the
     * container itself, `geo::PlaneDataContainer`. Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto hitsPerPlane
     *   = geom->makePlaneData<std::vector<recob::Hit const*>>();
     *
     * for (recob::Hit const& hit: hits) {
     *   if (hit.WireID()) hitsPerPlane[hit.WireID()].push_back(&hit);
     * } // for
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * where the container will be filled with pointers to all hits on the given
     * wire plane (wire IDs are implicitly converted into plane IDs in the index
     * `operator[]` call).
     */
    template <typename T>
    PlaneDataContainer<T> makePlaneData() const
    {
      return {Ncryostats(), MaxTPCs(), MaxPlanes()};
    }

    /**
     * @brief Returns a container with one entry per wire plane.
     * @tparam T type of data in the container
     * @param defValue the initial value of all elements in the container
     * @return a container with one default-constructed `T` per plane
     * @see `geo::PlaneDataContainer`
     *
     * This function operates as `makePlaneData() const`, except that copies
     * the specified value into all the entries of the container. Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto nHitsPerPlane = geom->makePlaneData(0U);
     *
     * for (recob::Hit const& hit: hits) {
     *   if (hit.WireID()) ++(hitsPerPlane[hit.WireID()]);
     * } // for
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template <typename T>
    PlaneDataContainer<T> makePlaneData(T const& defValue) const
    {
      return {Ncryostats(), MaxTPCs(), MaxPlanes(), defValue};
    }

    //@{
    /**
     * @brief Returns the total number of planes in the specified TPC
     * @param tpcid TPC ID
     * @return number of planes in specified TPC, or 0 if no TPC found
     *
     * The NElements() and NSiblingElements() methods are overloaded and their
     * return depends on the type of ID.
     *
     * @todo Change return type to size_t
     */
    unsigned int Nplanes(TPCID const& tpcid) const
    {
      TPCGeo const* pTPC = GetElementPtr(tpcid);
      return pTPC ? pTPC->NElements() : 0;
    }
    unsigned int NElements(TPCID const& tpcid) const { return Nplanes(tpcid); }
    unsigned int NSiblingElements(PlaneID const& planeid) const { return Nplanes(planeid); }
    //@}

    /**
     * @brief Returns the number of views (different wire orientations)
     *
     * Returns the number of different views, or wire orientations, in the
     * detector.
     *
     * The function assumes that all TPCs in all cryostats of a detector have
     * the same number of planes, which should be a safe assumption.
     *
     * @todo Change return type to size_t
     */
    unsigned int Nviews() const;

    //
    // access
    //

    //@{
    /**
     * @brief Returns whether we have the specified plane
     *
     * The HasElement() method is overloaded and its meaning depends on the type
     * of ID.
     *
     */
    bool HasPlane(PlaneID const& planeid) const
    {
      TPCGeo const* pTPC = TPCPtr(planeid);
      return pTPC ? pTPC->HasPlane(planeid) : false;
    }
    bool HasElement(PlaneID const& planeid) const { return HasPlane(planeid); }
    //@}

    //@{
    /**
     * @brief Returns the specified wire
     * @param planeid ID of the plane
     * @param p plane number within the TPC
     * @param tpc TPC number within the cryostat
     * @param cstat number of cryostat
     * @return a constant reference to the specified plane
     * @throw cet::exception (`GeometryCore` category) if cryostat not present
     * @throw cet::exception (`TPCOutOfRange` category) if no such TPC
     * @throw cet::exception (`PlaneOutOfRange` category) if no such plane
     *
     * The GetElement() method is overloaded and its return depends on the type
     * of ID.
     *
     * @todo remove the version with integers
     */
    PlaneGeo const& Plane(unsigned int const p,
                          unsigned int const tpc = 0,
                          unsigned int const cstat = 0) const
    {
      return Plane(PlaneID(cstat, tpc, p));
    }
    PlaneGeo const& Plane(PlaneID const& planeid) const { return TPC(planeid).Plane(planeid); }
    PlaneGeo const& GetElement(PlaneID const& planeid) const { return Plane(planeid); }
    //@}

    //@{
    /**
     * @brief Returns the specified plane
     * @param planeid plane ID
     * @return a constant pointer to the specified plane, or nullptr if none
     *
     * The GetElementPtr() method is overloaded and its return depends on the
     * type of ID.
     */
    PlaneGeo const* PlanePtr(PlaneID const& planeid) const
    {
      TPCGeo const* pTPC = TPCPtr(planeid);
      return pTPC ? pTPC->PlanePtr(planeid) : nullptr;
    }
    PlaneGeo const* GetElementPtr(PlaneID const& planeid) const { return PlanePtr(planeid); }
    //@}

    //
    // iterators
    //

    /// Initializes the specified ID with the ID of the first plane.
    void GetBeginID(PlaneID& id) const
    {
      GetBeginID(id.asTPCID());
      id.Plane = 0;
    }

    /// Initializes the specified ID with the invalid ID after the last plane.
    void GetEndID(PlaneID& id) const;

    /// Sets the ID to the ID after the specified one.
    /// @return whether the ID is actually valid (validity flag is also set)
    bool IncrementID(PlaneID& id) const; // inline implementation

    /// Returns the ID of the first plane of the specified cryostat.
    PlaneID GetBeginPlaneID(CryostatID const& id) const { return {GetBeginTPCID(id), 0}; }

    /// Returns the (possibly invalid) ID after the last plane of the specified
    /// cryostat.
    PlaneID GetEndPlaneID(CryostatID const& id) const;

    /// Returns the ID of the first plane of the specified TPC.
    PlaneID GetBeginPlaneID(TPCID const& id) const { return {id, 0}; }

    /// Returns the (possibly invalid) ID after the last plane of the specified
    /// TPC.
    PlaneID GetEndPlaneID(TPCID const& id) const;

    //
    // single object features
    //

    //@{
    /**
     * @brief Returns the distance between two planes
     * @param p1 index of the first plane
     * @param p2 index of the second plane
     * @param tpc tpc number within the cryostat
     * @param cstat cryostat number
     * @return distance between the planes
     *
     * @todo add a version with plane IDs
     * @todo deprecate this function
     * @todo add a default version for a given TPCID
     * @todo add a version with two plane indices for a given TPCID
     * @todo return the absolute value of the distance (makes the order unimportant)
     * @todo document what will happen (in the future methods) with planes on different TPCs
     */
    Length_t PlanePitch(TPCID const& tpcid,
                        PlaneID::PlaneID_t p1 = 0,
                        PlaneID::PlaneID_t p2 = 1) const;
    Length_t PlanePitch(PlaneID const& pid1, PlaneID const& pid2) const;
    Length_t PlanePitch(unsigned int p1 = 0,
                        unsigned int p2 = 1,
                        unsigned int tpc = 0,
                        unsigned int cstat = 0) const;
    //@}

    /**
     * @brief Returns the view (wire orientation) on the channels of specified TPC plane
     * @param plane TPC plane ID
     * @return the type of signal on the specified plane, or geo::kUnknown
     */
    View_t View(PlaneID const& pid) const;

    /**
     * @brief Returns the type of signal on the channels of specified TPC plane
     * @param plane TPC plane ID
     * @return the type of signal on the specified plane, or geo::kMysteryType
     *
     * Assumes that all the channels on the plane have the same signal type.
     *
     * @todo verify that kMysteryType is returned on invalid plane
     */
    SigType_t SignalType(PlaneID const& pid) const;

    /// @} Plane access and information

    /// @name Wire access and information
    /// @{

    //
    // group features
    //

    /**
     * @brief Returns the total number of wires in the specified plane
     * @param p plane number within the TPC
     * @param tpc tpc number within the cryostat
     * @param cstat cryostat number
     *
     * @todo Make all the arguments mandatory (as PlaneID)
     * @todo Change return type to size_t
     * @todo what happens if it does not exist?
     */
    unsigned int Nwires(unsigned int p, unsigned int tpc = 0, unsigned int cstat = 0) const
    {
      return Nwires(PlaneID(cstat, tpc, p));
    }

    //@{
    /**
     * @brief Returns the total number of wires in the specified plane
     * @param planeid plane ID
     * @return number of wires in specified plane, or 0 if no plane found
     *
     * The NElements() and NSiblingElements() methods are overloaded and their
     * return depends on the type of ID.
     *
     * @todo Change return type to size_t
     */
    unsigned int Nwires(PlaneID const& planeid) const
    {
      PlaneGeo const* pPlane = GetElementPtr(planeid);
      return pPlane ? pPlane->NElements() : 0;
    }
    unsigned int NElements(PlaneID const& planeid) const { return Nwires(planeid); }
    unsigned int NSiblingElements(WireID const& wireid) const { return Nwires(wireid); }

    /// Returns the largest number of wires among all planes in this detector
    unsigned int MaxWires() const;

    //@}

    //
    // access
    //

    //@{
    /**
     * @brief Returns whether we have the specified wire
     *
     * The HasElement() method is overloaded and its meaning depends on the type
     * of ID.
     */
    bool HasWire(WireID const& wireid) const
    {
      PlaneGeo const* pPlane = PlanePtr(wireid);
      return pPlane ? pPlane->HasWire(wireid) : false;
    }
    bool HasElement(WireID const& wireid) const { return HasWire(wireid); }
    //@}

    //@{
    /**
     * @brief Returns the specified wire
     * @param wireid wire ID
     * @return a constant pointer to the specified wire, or nullptr if none
     *
     * The GetElementPtr() method is overloaded and its return depends on the
     * type of ID.
     */
    WireGeo const* WirePtr(WireID const& wireid) const
    {
      PlaneGeo const* pPlane = PlanePtr(wireid);
      return pPlane ? pPlane->WirePtr(wireid) : nullptr;
    } // WirePtr()
    WireGeo const* GetElementPtr(WireID const& wireid) const { return WirePtr(wireid); }
    //@}

    //@{
    /**
     * @brief Returns the specified wire
     * @param wireid ID of the wire
     * @return a constant reference to the specified wire
     * @throw cet::exception if not found
     *
     * The GetElement() method is overloaded and its return depends on the type
     * of ID.
     */
    WireGeo const& Wire(WireID const& wireid) const { return Plane(wireid).Wire(wireid); }
    WireGeo const& WireIDToWireGeo(WireID const& wireid) const { return Wire(wireid); }
    WireGeo const& GetElement(WireID const& wireid) const { return Wire(wireid); }
    //@}

    //
    // iterators
    //

    /// Initializes the specified ID with the ID of the first wire.
    void GetBeginID(WireID& id) const
    {
      GetBeginID(id.asPlaneID());
      id.Wire = 0;
    }

    /// Initializes the specified ID with the invalid ID after the last wire.
    void GetEndID(WireID& id) const;

    /// Sets the ID to the ID after the specified one.
    /// @return whether the ID is actually valid (validity flag is also set)
    bool IncrementID(WireID& id) const; // inline implementation

    /// Returns the ID of the first wire in the specified cryostat.
    WireID GetBeginWireID(CryostatID const& id) const { return {GetBeginPlaneID(id), 0}; }

    /// Returns the (possibly invalid) ID after the last wire in the specified
    /// cryostat.
    WireID GetEndWireID(CryostatID const& id) const;

    /// Returns the ID of the first wire of the specified TPC.
    WireID GetBeginWireID(TPCID const& id) const { return {PlaneID(id, 0), 0}; }

    /// Returns the (possibly invalid) ID after the last wire of the specified
    /// TPC.
    WireID GetEndWireID(TPCID const& id) const;

    /// Returns the ID of the first wire of the specified wire plane.
    WireID GetBeginWireID(PlaneID const& id) const { return {id, 0}; }

    /// Returns the (possibly invalid) ID after the last wire of the specified
    /// wire plane.
    WireID GetEndWireID(PlaneID const& id) const;

    //
    // single object features
    //

    //@{
    /**
     * @brief Returns the distance between two consecutive wires.
     * @param p plane number within the TPC
     * @param tpc tpc number within the cryostat
     * @param cstat cryostat number
     * @return the distance between the two wires
     *
     * @note The current geometry assumptions imply that wire pitch is constant
     *       between all wires on the same wire plane. This is an assumption
     *       non-trivial to remove.
     *
     * @todo add a version with wire IDs
     * @todo deprecate this function
     * @todo document what will happen (in the future methods) with wires on different planes
     *
     */
    Length_t WirePitch(PlaneID const& planeid) const;
    Length_t WirePitch(unsigned int plane = 0, unsigned int tpc = 0, unsigned int cstat = 0) const
    {
      return WirePitch(PlaneID(cstat, tpc, plane));
    }
    //@}

    /**
     * @brief Returns the distance between two wires in the specified view
     * @param w1 index of the first wire
     * @param w2 index of the second wire
     * @param p plane number within the TPC
     * @param tpc tpc number within the cryostat
     * @param cstat cryostat number
     * @return the distance between the two wires
     *
     * This method assumes that all the wires on all the planes on the specified
     * view of all TPCs have the same pitch.
     */
    Length_t WirePitch(View_t view) const;

    //@{
    /**
     * @brief Returns the angle of the wires in the specified view from vertical
     * @param view the view
     * @param TPC the index of the TPC in the specified cryostat
     * @param Cryo the cryostat
     * @param tpcid ID of the TPC
     * @return the angle [radians]
     * @throw cet::exception ("GeometryCore" category) if no such view
     *
     * The angle is defined as in WireGeo::ThetaZ().
     *
     * This method assumes all wires in the view have the same angle (it queries
     * for the first).
     *
     * @deprecated This does not feel APA-ready
     */
    double WireAngleToVertical(View_t view, TPCID const& tpcid) const;
    double WireAngleToVertical(View_t view, int TPC = 0, int Cryo = 0) const
    {
      return WireAngleToVertical(view, TPCID(Cryo, TPC));
    }
    //@}

    /// @} Wire access and information

    /**
     * @name Wire geometry queries
     *
     * Please note the differences between functions:
     * ChannelsIntersect(), WireIDsIntersect() and IntersectionPoint()
     * all calculate wires intersection using the same equation.
     * ChannelsIntersect() and WireIdsIntersect() will return true
     * if the two wires cross, return false if they don't.
     * IntersectionPoint() does not check if the two wires cross.
     */
    /// @{

    //
    // simple geometry queries
    //

    /**
     * @brief Fills two arrays with the coordinates of the wire end points
     * @param wireid ID of the wire
     * @param xyzStart (output) an array with the start coordinate
     * @param xyzEnd (output) an array with the end coordinate
     * @throws cet::exception wire not present
     *
     * The starting point is the wire end with lower z coordinate.
     *
     * @deprecated use the wire ID interface instead (but note that it does not
     *             sort the ends)
     */
    void WireEndPoints(WireID const& wireid, double* xyzStart, double* xyzEnd) const;

    /**
     * @brief Fills two arrays with the coordinates of the wire end points
     * @param cstat cryostat number
     * @param tpc tpc number within the cryostat
     * @param plane plane number within the TPC
     * @param wire wire number within the plane
     * @param xyzStart (output) an array with the start coordinate
     * @param xyzEnd (output) an array with the end coordinate
     * @throws cet::exception wire not present
     *
     * The starting point is the wire end with lower z coordinate.
     *
     * @deprecated use the wire ID interface instead (but note that it does not
     *             sort the ends)
     */
    void WireEndPoints(unsigned int cstat,
                       unsigned int tpc,
                       unsigned int plane,
                       unsigned int wire,
                       double* xyzStart,
                       double* xyzEnd) const
    {
      WireEndPoints(WireID(cstat, tpc, plane, wire), xyzStart, xyzEnd);
    }

    //@{
    /**
     * @brief Returns a segment whose ends are the wire end points
     * @param wireid ID of the wire
     * @return a segment whose ends are the wire end points
     * @throws cet::exception wire not present
     *
     * The start and end are assigned as returned from the geo::WireGeo object.
     * The rules for this assignment are documented in that class.
     *
     * @deprecated use the wire ID interface instead (but note that it does not
     *             sort the ends)
     */
    template <typename Point>
    Segment<Point> WireEndPoints(WireID const& wireID) const;
    Segment<DefaultPoint_t> WireEndPoints(WireID const& wireID) const
    {
      return WireEndPoints<DefaultPoint_t>(wireID);
    }

    //@}

    //
    // closest wire
    //

    /**
     * @brief Returns the ID of wire closest to position in the specified TPC.
     * @param point the point to be tested [cm]
     * @param planeid ID of the plane
     * @return the ID of the wire, or an invalid wire ID
     * @see `geo::PlaneGeo::ClosestWireID()`
     * @bug Instead of returning an invalid wire ID, an exception is thrown!
     *
     * If the nearest wire is not closer than half a wire pitch, the result is
     * marked invalid. The returned (invalid) ID will contain the non-existing
     * wire that would be the nearest, if it existed.
     *
     * If the wire ID is invalid and the existing closest wire is desired,
     * a possible solution is (when the BUG will be solved):
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::WireID wireID = geom->NearestWireID(point, planeID);
     * if (!wireID) wireID = geom->Plane(planeID).ClosestWireID(wireID);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Note however that this will execute plane lookup twice, and a more
     * efficient approach would be to ask the plane everything directly:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::PlaneGeo const& plane = geom->Plane(planeID);
     * geo::WireID wireID = plane.NearestWireID(point);
     * if (!wireID) wireID = plane.ClosestWireID(wireID);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     * Until the BUG is fixed, the actual working code is:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * geo::PlaneGeo const& plane = geom->Plane(planeID);
     * geo::WireID wireID;
     * try {
     *   wireID = plane.NearestWireID(point);
     * }
     * catch (geo::InvalidWireError const& e) {
     *   if (!e.hasSuggestedWire()) throw;
     *   wireID = plane.ClosestWireID(e.suggestedWireID());
     * }
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    WireID NearestWireID(Point_t const& point, PlaneID const& planeid) const;

    //@{
    /**
     * @brief Returns the index of the nearest wire to the specified position
     * @param pos world coordinates of the position (it will be projected)
     * @param planeid ID of the plane
     * @return an index interpolation between the two nearest wires
     * @see ChannelMapAlg::WireCoordinate()
     *
     * Respect to NearestWireID(), this method returns a real number,
     * representing a continuous coordinate in the wire axis, with the round
     * values corresponding to the actual wires.
     */
    Length_t WireCoordinate(Point_t const& pos, PlaneID const& planeid) const;
    //@}

    //
    // wire intersections
    //

    bool IntersectLines(double A_start_x,
                        double A_start_y,
                        double A_end_x,
                        double A_end_y,
                        double B_start_x,
                        double B_start_y,
                        double B_end_x,
                        double B_end_y,
                        double& x,
                        double& y) const
    {
      return unrelated::IntersectLines(
        A_start_x, A_start_y, A_end_x, A_end_y, B_start_x, B_start_y, B_end_x, B_end_y, x, y);
    }

    bool IntersectSegments(double A_start_x,
                           double A_start_y,
                           double A_end_x,
                           double A_end_y,
                           double B_start_x,
                           double B_start_y,
                           double B_end_x,
                           double B_end_y,
                           double& x,
                           double& y) const
    {
      return unrelated::IntersectSegments(
        A_start_x, A_start_y, A_end_x, A_end_y, B_start_x, B_start_y, B_end_x, B_end_y, x, y);
    }

    //@{
    /**
     * @brief Computes the intersection between two wires.
     * @param wid1 ID of the first wire
     * @param wid2 ID of the other wire
     * @param[out] intersection the intersection point (global coordinates)
     * @return whether an intersection was found inside the TPC the wires belong
     * @see `geo::WiresIntersection()`, `geo::LineClosestPoint()`
     *
     * The wires identified by `wid1` and `wid2` are intersected, and the
     * 3D intersection point is written into the `intersection` parameter.
     * The "intersection" point is actually the point belonging to the first
     * wire (`wid2`) which is the closest (in Euclidean 3D metric) to the second
     * wire.
     *
     * The intersection is computed only if the wires belong to different planes
     * of the same TPC. If that is not the case (i.e. they belong to different
     * TPC or cryostat, or if they belong to the same plane), `false` is
     * returned and `intersection` is set with all components to infinity
     * (`std::numeric_limits<>::infinity()`).
     *
     * When the intersection is computed, it is always stored in the
     * `intersection` output parameter. Return value is `true` if this
     * intersection lies within the physical boundaries first wire, while it is
     * instead `false` if it lies on the extrapolation of the wire direction,
     * but not within the wire physical extension.
     *
     * To test that the result is not infinity (nor NaN), use
     * `geo::vect::isfinite(intersection)` etc.
     *
     * @note If `geo::WireGeo` objects are already available, using instead
     *       the free function `geo::WiresIntersection()` or the method
     *       `geo::WireGeo::IntersectionWith()` is faster (and _recommended_).
     *       For purely geometric intersection, `geo::LineClosestPoint()` is
     *       also available.
     */
    bool WireIDsIntersect(WireID const& wid1, WireID const& wid2, Point_t& intersection) const;
    bool WireIDsIntersect(WireID const& wid1, WireID const& wid2, TVector3& intersection) const;
    //@}

    //@{
    /**
     * @brief Computes the intersection between two wires.
     * @param wid1 ID of the first wire
     * @param wid2 ID of the other wire
     * @param widIntersect (output) the coordinate of the intersection point
     * @return whether an intersection was found within the TPC
     *
     * The "intersection" refers to the projection of the wires into the same
     * @f$ x = 0 @f$ plane.
     * Wires are assumed to have at most one intersection.
     * If wires are parallel, `widIntersect` will have the two components set to
     * infinity (`std::numeric_limits<>::infinity()`) and the TPC number set to
     * invalid (`geo::TPCID::InvalidID`). Also, `false` is returned.
     * If the intersection is outside the TPC, `false` is also returned, but the
     * `widIntersect` will contain the coordinates of that intersection. The TPC
     * number is still set to invalid, although the intersection _might_ belong
     * to a valid TPC somewhere else.
     *
     *
     * @deprecated This method uses arbitrary assumptions and should not be
     *             used. Use the interface returning a full vector instead.
     */
    bool WireIDsIntersect(WireID const& wid1,
                          WireID const& wid2,
                          WireIDIntersection& widIntersect) const;
    //@}

    /**
     * @brief Returns the intersection point of two wires
     * @param wid1 ID of the first wire
     * @param wid2 ID of the other wire
     * @param y (output) y coordinate of the intersection point
     * @param z (output) z coordinate of the intersection point
     * @return whether an intersection was found within the TPC
     * @see WireIDsIntersect()
     *
     * The behaviour of this method reflects the one of `WireIDsIntersect()`,
     * which supersedes this one.
     *
     * To test if the result is infinity, use e.g. `std::isfinite(y)`.
     *
     * @deprecated This method uses arbitrary assumptions and should not be
     *             used. Use `WireIDsIntersect()` returning a vector, instead.
     */
    bool IntersectionPoint(WireID const& wid1, WireID const& wid2, double& y, double& z) const;

    /**
     * @brief Returns the intersection point of two wires
     * @param wire1 wire index of the first wire
     * @param wire2 wire index of the other wire
     * @param plane1 plane index of the first wire
     * @param plane2 plane index of the other wire
     * @param cstat cryostat number
     * @param tpc tpc number within the cryostat where the planes belong
     * @param y (output) y coordinate of the intersection point
     * @param z (output) z coordinate of the intersection point
     * @return whether an intersection was found
     *
     * No check is performed, not any information provided, about the validity
     * of the result.
     *
     * @deprecated This method uses arbitrary assumptions and should not be
     *             used. Use `WireIDsIntersect()` returning a vector, instead.
     */
    bool IntersectionPoint(unsigned int wire1,
                           unsigned int wire2,
                           unsigned int plane1,
                           unsigned int plane2,
                           unsigned int cstat,
                           unsigned int tpc,
                           double& y,
                           double& z) const
    {
      return IntersectionPoint(
        WireID(cstat, tpc, plane1, wire1), WireID(cstat, tpc, plane2, wire2), y, z);
    }

    /**
     * @brief Returns the plane that is not in the specified arguments
     * @param pid1 a plane
     * @param pid2 another plane
     * @return the ID to the third plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     * @throws cet::exception (category: "GeometryCore") if pid1 and pid2 match
     *
     * This function requires a geometry with exactly three planes.
     * If the two input planes are not on the same TPC, the result is undefined.
     */
    PlaneID ThirdPlane(PlaneID const& pid1, PlaneID const& pid2) const;

    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param pid1 ID of the plane of the first slope
     * @param slope1 slope as seen on the first plane
     * @param pid2 ID of the plane of the second slope
     * @param slope2 slope as seen on the second plane
     * @param output_plane ID of the plane on which to calculate the slope
     * @return the slope on the third plane, or -999. if slope would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if input planes match
     *
     * Given a slope as projected in two planes, returns the slope as projected
     * in the specified output plane.
     * The slopes are defined in uniform units; they should be computed as
     * distance ratios (or tangent of a geometrical angle; the formula is still
     * valid using dt/dw directly in case of equal wire pitch in all planes
     * and uniform drift velocity.
     */
    double ThirdPlaneSlope(PlaneID const& pid1,
                           double slope1,
                           PlaneID const& pid2,
                           double slope2,
                           PlaneID const& output_plane) const;

    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param pid1 ID of the plane of the first slope
     * @param slope1 slope as seen on the first plane
     * @param pid2 ID of the plane of the second slope
     * @param slope2 slope as seen on the second plane
     * @return the slope on the third plane, or -999. if slope would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if same plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     *
     * Given a slope as projected in two planes, returns the slope as projected
     * in the third plane.
     * This function is a shortcut assuming exactly three wire planes in the
     * TPC, in which case the output plane is chosen as the one that is neither
     * of the input planes.
     */
    double ThirdPlaneSlope(PlaneID const& pid1,
                           double slope1,
                           PlaneID const& pid2,
                           double slope2) const;

    //@{
    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param plane1 index of the plane of the first slope
     * @param slope1 slope as seen on the first plane
     * @param plane2 index of the plane of the second slope
     * @param slope2 slope as seen on the second plane
     * @param tpcid TPC where the two planes belong
     * @return the slope on the third plane, or -999. if slope would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if same plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     *
     * Given a slope as projected in two planes, returns the slope as projected
     * in the third plane.
     */
    double ThirdPlaneSlope(PlaneID::PlaneID_t plane1,
                           double slope1,
                           PlaneID::PlaneID_t plane2,
                           double slope2,
                           TPCID const& tpcid) const
    {
      return ThirdPlaneSlope(PlaneID(tpcid, plane1), slope1, PlaneID(tpcid, plane2), slope2);
    }
    double ThirdPlaneSlope(unsigned int plane1,
                           double slope1,
                           unsigned int plane2,
                           double slope2,
                           unsigned int tpc,
                           unsigned int cstat) const
    {
      return ThirdPlaneSlope(plane1, slope1, plane2, slope2, TPCID(cstat, tpc));
    }
    //@}

    /**
     * @brief Returns dT/dW on the third plane, given it in the other two
     * @param pid1 ID of the plane of the first dT/dW
     * @param dTdW1 dT/dW as seen on the first plane
     * @param pid2 ID of the plane of the second dT/dW
     * @param dTdW2 dT/dW  as seen on the second plane
     * @param output_plane ID of the plane on which to calculate the slope
     * @return dT/dW on the third plane, or -999. if dT/dW would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if same plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     *
     * Given a dT/dW as projected in two planes, returns the dT/dW as projected
     * in the third plane.
     * The dT/dW are defined in time ticks/wide number units.
     */
    double ThirdPlane_dTdW(PlaneID const& pid1,
                           double slope1,
                           PlaneID const& pid2,
                           double slope2,
                           PlaneID const& output_plane) const;

    /**
     * @brief Returns dT/dW on the third plane, given it in the other two
     * @param pid1 ID of the plane of the first dT/dW
     * @param dTdW1 dT/dW as seen on the first plane
     * @param pid2 ID of the plane of the second dT/dW
     * @param dTdW2 dT/dW  as seen on the second plane
     * @return dT/dW on the third plane, or -999. if dT/dW would be infinity
     * @throws cet::exception (category: "GeometryCore") if different TPC
     * @throws cet::exception (category: "GeometryCore") if same plane
     * @throws cet::exception (category: "GeometryCore") if other than 3 planes
     *
     * Given a dT/dW as projected in two planes, returns the dT/dW as projected
     * in the third plane.
     * This function is a shortcut assuming exactly three wire planes in the
     * TPC, in which case the output plane is chosen as the one that is neither
     * of the input planes.
     */
    double ThirdPlane_dTdW(PlaneID const& pid1,
                           double slope1,
                           PlaneID const& pid2,
                           double slope2) const;

    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param angle1 angle or the wires on the first plane
     * @param slope1 slope as observed on the first plane
     * @param angle2 angle or the wires on the second plane
     * @param slope2 slope as observed on the second plane
     * @param angle_target angle or the wires on the target plane
     * @return the slope as measure on the third plane, or 999 if infinity
     *
     * This function will return a small slope if both input slopes are small.
     */
    static double ComputeThirdPlaneSlope(double angle1,
                                         double slope1,
                                         double angle2,
                                         double slope2,
                                         double angle_target);

    /**
     * @brief Returns the slope on the third plane, given it in the other two
     * @param angle1 angle or the wires on the first plane
     * @param pitch1 wire pitch on the first plane
     * @param dTdW1 slope in dt/dw units as observed on the first plane
     * @param angle2 angle or the wires on the second plane
     * @param pitch2 wire pitch on the second plane
     * @param dTdW2 slope in dt/dw units as observed on the second plane
     * @param angle_target angle or the wires on the target plane
     * @param pitch_target wire pitch on the target plane
     * @return dt/dw slope as measured on the third plane, or 999 if infinity
     *
     * The input slope must be specified in dt/dw non-homogeneous coordinates.
     *
     * This function will return a small slope if both input slopes are small.
     */
    static double ComputeThirdPlane_dTdW(double angle1,
                                         double pitch1,
                                         double dTdW1,
                                         double angle2,
                                         double pitch2,
                                         double dTdW2,
                                         double angle_target,
                                         double pitch_target);

    /// @} Wire geometry queries

    /**
     * @name Optical detector geometry access and information
     * @anchor GeometryCoreOpDetGeometry
     * @see @ref GeometryCoreOpDetChannel "optical detector channel information"
     *
     * There are a number of ways to identify an optical detector or channel:
     *
     * * geometric:
     *     * cryostat (e.g. `geo::CryostatID`) and relative optical detector
     *       number within it
     *     * unique optical detector number
     * * readout:
     *     * optical detector channel
     *     * "hardware" channel
     *
     * And they all should be better documented!
     */
    /// @{

    //
    // group features
    //

    /// Number of OpDets in the whole detector
    unsigned int NOpDets() const;

    //
    // access
    //
    /**
     * @brief Returns the `geo::OpDetGeo` object for the given channel number.
     * @param OpChannel optical detector unique channel number
     * @see GeometryCoreOpDetGeometry "optical detector identification"
     */
    OpDetGeo const& OpDetGeoFromOpChannel(unsigned int OpChannel) const;

    /**
     * @brief Returns the `geo::OpDetGeo` object for the given detector number.
     * @param OpDet optical detector unique number
     * @see GeometryCoreOpDetGeometry "optical detector identification"
     */
    OpDetGeo const& OpDetGeoFromOpDet(unsigned int OpDet) const;

    //@{
    /**
     * @brief Find the nearest OpChannel to some point
     * @param xyz point to be queried, in world coordinates
     * @return the nearest OpChannel to the point,
     *         or `std::numeric_limits<unsigned int>::max()` if invalid point
     *
     * @deprecated This method does not tell in which cryostat the detector is;
     *             use `geo::CryostatGeo::GetClosestOpDet()` instead
     *             (find the cryostat with `PositionToCryostatPtr()`).
     *
     */
    unsigned int GetClosestOpDet(Point_t const& point) const;
    unsigned int GetClosestOpDet(double const* point) const;
    //@}

    //
    // object description
    //

    /**
     * @brief Returns gdml string which gives sensitive opdet name
     * @param c ID of the cryostat the detector is in
     *
     * This name is defined in the geometry (GDML) description.
          *
     * @todo Change to use CryostatID
     */
    std::string OpDetGeoName(unsigned int c = 0) const;

    /// @} Optical detector access and information

    /// @name Auxiliary detectors access and information
    /// @{

    /// @todo use a AutDetID_t instead of unsigned int?

    //
    // group features
    //

    /**
     * @brief Returns the number of auxiliary detectors
     *
     * This method returns the total number of scintillator paddles
     * (Auxiliary Detectors aka AuxDet) outside of the cryostat
     *
     * @todo Change return type to size_t
     */
    unsigned int NAuxDets() const { return AuxDets().size(); }

    /**
     * @brief Returns the number of sensitive components of auxiliary detector
     * @param aid ID of the auxiliary detector
     * @return number of sensitive components in the auxiliary detector aid
     * @thrws cet::exception (category "Geometry") if aid does not exist
     */
    unsigned int NAuxDetSensitive(size_t const& aid) const;

    //
    // access
    //

    /**
     * @brief Returns the specified auxiliary detector
     * @param ad the auxiliary detector index
     * @return a constant reference to the specified auxiliary detector
     *
     * @todo what happens if it does not exist?
     * @todo remove the default parameter?
     */
    AuxDetGeo const& AuxDet(unsigned int const ad = 0) const;

    /**
     * @brief Returns the index of the auxiliary detector at specified location.
     * @param point location to be tested
     * @param tolerance tolerance (cm) for matches. Default 0
     * @return the index of the detector, or
     *        `std::numeric_limits<unsigned int>::max()` if no detector is there
     *
     * @bug Actually, an exception is thrown.
     */
    unsigned int FindAuxDetAtPosition(Point_t const& point, double tolerance = 0) const;

    /**
     * @brief Fills the indices of the sensitive auxiliary detector at location
     * @param point location to be tested
     * @param adg _(output)_ auxiliary detector index
     * @param sv _(output)_ sensitive volume index
     * @param tolerance tolerance (cm) for matches. Default 0.
     */
    void FindAuxDetSensitiveAtPosition(Point_t const& point,
                                       std::size_t& adg,
                                       std::size_t& sv,
                                       double tolerance = 0) const;

    /**
     * @brief Returns the auxiliary detector at specified location
     * @param point location to be tested
     * @param ad _(output)_ the auxiliary detector index
     * @param tolerance tolerance (cm) for matches. Default 0.
     * @return constant reference to AuxDetGeo object of the auxiliary detector
     *
     * @todo what happens if it does not exist?
     */
    AuxDetGeo const& PositionToAuxDet(Point_t const& point,
                                      unsigned int& ad,
                                      double tolerance = 0) const;

    /**
     * @brief Returns the auxiliary detector at specified location
     * @param point location to be tested
     * @param ad _(output)_ the auxiliary detector index
     * @param sv _(output)_ the auxiliary detector sensitive volume index
     * @param tolerance tolerance (cm) for matches. Default 0.
     * @return reference to AuxDetSensitiveGeo object of the auxiliary detector
     *
     * @todo what happens if it does not exist?
     */
    const AuxDetSensitiveGeo& PositionToAuxDetSensitive(Point_t const& point,
                                                        size_t& ad,
                                                        size_t& sv,
                                                        double tolerance = 0) const;

    const AuxDetGeo& ChannelToAuxDet(
      std::string const& auxDetName,
      uint32_t const& channel) const; // return the AuxDetGeo for the given detector
                                      // name and channel

    const AuxDetSensitiveGeo& ChannelToAuxDetSensitive(
      std::string const& auxDetName,
      uint32_t const& channel) const; // return the AuxDetSensitiveGeo for the given

    /// @} Auxiliary detectors access and information

    /// @name TPC readout channels and views
    /// @{

    //
    // group features
    //

    /// Returns the number of TPC readout channels in the detector
    unsigned int Nchannels() const;

    /// @brief Returns the number of channels in the specified ROP
    /// @return number of channels in the specified ROP, 0 if non-existent
    unsigned int Nchannels(readout::ROPID const& ropid) const;

    /// @brief Returns an std::vector<ChannelID_t> in all TPCs in a TPCSet
    std::vector<raw::ChannelID_t> ChannelsInTPCs() const;
    //
    /**
     * @brief Returns a list of possible views in the detector.
     * @return the set of views
     */
    std::set<View_t> const& Views() const { return allViews; }

    //
    // access
    //

    /**
     * @brief Returns whether the specified channel exists and is valid
     * @param channel the ID of the channel
     * @return whether the specified channel exists
     *
     * A channel is defined as existing and valid if its ID is not invalid and
     * if the channel is physical.
     */
    bool HasChannel(raw::ChannelID_t channel) const;

    //@{
    /**
     * @brief Returns the ID of the TPC channel connected to the specified wire
     * @param plane the number of plane
     * @param wire the number of wire
     * @param tpc the number of TPC
     * @param cryostat the number of cryostat
     * @param wireid the ID of the wire
     * @return the ID of the channel, or raw::InvalidChannelID if invalid wire
     *
     * @todo Verify the raw::InvalidChannelID part
     * @todo remove the integers version
     */
    raw::ChannelID_t PlaneWireToChannel(WireID const& wireid) const;
    raw::ChannelID_t PlaneWireToChannel(unsigned int const plane,
                                        unsigned int const wire,
                                        unsigned int const tpc = 0,
                                        unsigned int const cstat = 0) const
    {
      return PlaneWireToChannel(WireID(cstat, tpc, plane, wire));
    }
    //@}

    //
    // single object features
    //

    /**
     * @brief Returns the type of signal on the specified TPC channel
     * @param channel TPC channel ID
     * @return the type of signal on the specified channel, or geo::kMysteryType
     *
     * @todo verify that kMysteryType is returned on invalid channel
     */
    SigType_t SignalType(raw::ChannelID_t const channel) const;

    /**
     * @brief Returns the view (wire orientation) on the specified TPC channel
     * @param channel TPC channel ID
     * @return the type of signal on the specified channel, or geo::kUnknown
     *
     * The view of the readout plane `channel` belongs to is returned, as in
     * `View(readout::ROPID const&) const`.
     */
    View_t View(raw::ChannelID_t const channel) const;

    /**
     * @brief Returns a list of wires connected to the specified TPC channel
     * @param channel TPC channel ID
     * @return vector containing the ID of all the connected wires
     * @throws cet::exception (category: "Geometry") if non-existent channel
     */
    std::vector<WireID> ChannelToWire(raw::ChannelID_t const channel) const;

    /// Returns the ID of the ROP the channel belongs to
    /// @throws cet::exception (category: "Geometry") if non-existent channel
    readout::ROPID ChannelToROP(raw::ChannelID_t channel) const;

    //
    // geometry queries
    //

    /**
     * @brief Returns the ID of the channel nearest to the specified position
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param planeid ID of the wire plane the channel must belong to
     * @return the ID of the channel, or `raw::InvalidChannelID` if invalid wire
     * @bug on invalid wire, a `geo::InvalidWireError` exception is thrown
     *
     */
    raw::ChannelID_t NearestChannel(Point_t const& worldLoc, PlaneID const& planeid) const;

    //@{
    /**
     * @brief Returns the ID of the channel nearest to the specified position
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @param PlaneNo the number of plane
     * @param TPCNo the number of TPC
     * @param cstat the number of cryostat
     * @return the ID of the channel, or raw::InvalidChannelID if invalid wire
     * @bug on invalid wire, a `geo::InvalidWireError` exception is thrown
     *
     * The different versions allow different way to provide the position.
     *
     * @todo remove the integers version
     */
    raw::ChannelID_t NearestChannel(const double worldLoc[3], PlaneID const& planeid) const;
    raw::ChannelID_t NearestChannel(std::vector<double> const& worldLoc,
                                    PlaneID const& planeid) const;
    raw::ChannelID_t NearestChannel(const TVector3& worldLoc, PlaneID const& planeid) const
    {
      return NearestChannel(vect::toPoint(worldLoc), planeid);
    }
    raw::ChannelID_t NearestChannel(const double worldLoc[3],
                                    unsigned int const PlaneNo,
                                    unsigned int const TPCNo = 0,
                                    unsigned int const cstat = 0) const
    {
      return NearestChannel(worldLoc, PlaneID(cstat, TPCNo, PlaneNo));
    }
    raw::ChannelID_t NearestChannel(std::vector<double> const& worldLoc,
                                    unsigned int const PlaneNo,
                                    unsigned int const TPCNo = 0,
                                    unsigned int const cstat = 0) const
    {
      return NearestChannel(worldLoc, PlaneID(cstat, TPCNo, PlaneNo));
    }
    raw::ChannelID_t NearestChannel(const TVector3& worldLoc,
                                    unsigned int const PlaneNo,
                                    unsigned int const TPCNo = 0,
                                    unsigned int const cstat = 0) const
    {
      return NearestChannel(worldLoc, PlaneID(cstat, TPCNo, PlaneNo));
    }
    raw::ChannelID_t NearestChannel(Point_t const& worldLoc,
                                    unsigned int const PlaneNo,
                                    unsigned int const TPCNo = 0,
                                    unsigned int const cstat = 0) const
    {
      return NearestChannel(worldLoc, PlaneID(cstat, TPCNo, PlaneNo));
    }
    //@}

    /**
     * @brief Returns an intersection point of two channels
     * @param c1 one channel ID
     * @param c2 the other channel ID
     * @param y (output) y coordinate of the intersection
     * @param z (output) z coordinate of the intersection
     * @return whether a intersection point was found
     *
     * @todo what happens for channels from different TPCs?
     * @todo what happens for channels with multiple intersection points?
     *
     * @deprecated This is clearly not APA-aware
     */
    bool ChannelsIntersect(raw::ChannelID_t c1, raw::ChannelID_t c2, double& y, double& z) const;

    /// @} TPC readout channels

    /// @name TPC set information
    /// @{

    //
    // group features
    //

    //@{
    /**
     * @brief Returns the total number of TPC sets in the specified cryostat
     * @param cryoid cryostat ID
     * @return number of TPC sets in the cryostat, or 0 if no cryostat found
     *
     * The NSiblingElements() method is overloaded and its
     * return depends on the type of ID.
     */
    unsigned int NTPCsets(readout::CryostatID const& cryoid) const;
    unsigned int NSiblingElements(readout::TPCsetID const& tpcsetid) const
    {
      return NTPCsets(tpcsetid);
    }
    //@}

    /// Returns the largest number of TPC sets any cryostat in the detector has
    unsigned int MaxTPCsets() const;

    /**
     * @brief Returns a container with one entry per TPC set.
     * @tparam T type of data in the container
     * @return a container with one default-constructed `T` per TPC set
     * @see `readout::TPCsetDataContainer`
     *
     * The working assumption is that all cryostats have the same number of
     * TPC sets. It is always guaranteed that all existing TPC sets have an
     * entry in the container, although if the previous working assumption is
     * not satisfied there will be entries in the containers which are not
     * associated to a valid TPC set.
     *
     * The interface of the container is detailed in the documentation of the
     * container itself, `readout::TPCsetDataContainer`. Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto tracksPerTPCset
     *   = geom->makeTPCsetData<std::vector<recob::Track const*>>();
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template <typename T>
    readout::TPCsetDataContainer<T> makeTPCsetData() const
    {
      return {Ncryostats(), MaxTPCsets()};
    }

    /**
     * @brief Returns a container with one entry per TPC set.
     * @tparam T type of data in the container
     * @param defValue the initial value of all elements in the container
     * @return a container with a value `defValue` per each TPC set
     * @see `readout::TPCsetDataContainer`
     *
     * This function operates as `makeTPCsetData() const`, except that it copies
     * the specified value into all the entries of the container. Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto nTracksPerTPCset = geom->makeTPCsetData(0U);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template <typename T>
    readout::TPCsetDataContainer<T> makeTPCsetData(T const& defValue) const
    {
      return {Ncryostats(), MaxTPCsets(), defValue};
    }

    //
    // access
    //
    /// Returns whether we have the specified TPC set
    /// @return whether the TPC set is valid and exists
    bool HasTPCset(readout::TPCsetID const& tpcsetid) const;

    /// Returns whether we have the specified TPC set
    bool HasElement(readout::TPCsetID const& tpcsetid) const { return HasTPCset(tpcsetid); }

    /**
     * @brief Returns the ID of the TPC set at specified location
     * @param worldLoc 3D coordinates of the point (world reference frame)
     * @return the TPC set ID, or an invalid one if no TPC set is there
     */
    readout::TPCsetID FindTPCsetAtPosition(double const worldLoc[3]) const;

    //
    // mapping
    //
    /// Returns the ID of the TPC set tpcid belongs to
    readout::TPCsetID TPCtoTPCset(TPCID const& tpcid) const;

    /**
     * @brief Returns a list of ID of TPCs belonging to the specified TPC set
     * @param tpcsetid ID of the TPC set to convert into TPC IDs
     * @return the list of TPCs, empty if TPC set is invalid
     *
     * Note that the check is performed on the validity of the TPC set ID, that
     * does not necessarily imply that the TPC set specified by the ID actually
     * exists. Check the existence of the TPC set first (HasTPCset()).
     * Behaviour on valid, non-existent TPC set IDs is undefined.
     */
    std::vector<TPCID> TPCsetToTPCs(readout::TPCsetID const& tpcsetid) const;

    ///
    /// iterators
    ///

    /// Initializes the specified ID with the ID of the first TPC set
    void GetBeginID(readout::TPCsetID& id) const
    {
      GetBeginID(id.asCryostatID());
      id.TPCset = 0;
    }

    /// Initializes the specified ID with the invalid ID after the last TPC set
    void GetEndID(readout::TPCsetID& id) const
    {
      GetEndID(id.asCryostatID());
      id.TPCset = 0;
    }

    /// Sets the ID to the ID after the specified one.
    /// @return whether the ID is actually valid (validity flag is also set)
    bool IncrementID(readout::TPCsetID& id) const; // inline implementation

    /// Returns the ID of the first TPC set in the specified cryostat.
    readout::TPCsetID GetBeginTPCsetID(CryostatID const& id) const { return {id, 0}; }

    /// Returns the (possibly invalid) ID after the last TPC set of the
    /// specified cryostat.
    readout::TPCsetID GetEndTPCsetID(CryostatID const& id) const { return {id.Cryostat + 1, 0}; }

    /// @} TPC set information

    /// @name Readout plane information
    /// @{

    //
    // group features
    //

    //@{
    /**
     * @brief Returns the total number of ROP in the specified TPC set
     * @param tpcsetid TPC set ID
     * @return number of readout planes in the TPC set, or 0 if no TPC set found
     *
     * Note that this methods explicitly check the existence of the TPC set.
     *
     * The NSiblingElements() method is overloaded and its
     * return depends on the type of ID.
     */
    unsigned int NROPs(readout::TPCsetID const& tpcsetid) const;
    unsigned int NSiblingElements(readout::ROPID const& ropid) const { return NROPs(ropid); }
    //@}

    /// Returns the largest number of ROPs a TPC set in the detector has
    unsigned int MaxROPs() const;

    /**
     * @brief Returns a container with one entry per readout plane.
     * @tparam T type of data in the container
     * @return a container with one default-constructed `T` per readout plane
     * @see `readout::ROPDataContainer`
     *
     * The working assumption is that all cryostats have the same number of
     * TPC sets, and all TPC sets have the same number of readout planes.
     * It is always guaranteed that all existing readout planes have an entry
     * in the container, although if the previous working assumption is not
     * satisfied there will be entries in the container which are not
     * associated to a valid readout plane.
     *
     * The interface of the container is detailed in the documentation of the
     * container itself, `readout::ROPDataContainer`. Example of usage:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto hitsPerROP
     *   = geom->makeROPdata<std::vector<recob::Hit const*>>();
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template <typename T>
    readout::ROPDataContainer<T> makeROPdata() const
    {
      return {Ncryostats(), MaxTPCsets(), MaxROPs()};
    }

    /**
     * @brief Returns a container with one entry per readout plane.
     * @tparam T type of data in the container
     * @param defValue the initial value of all elements in the container
     * @return a container with one default-constructed `T` per readout plane
     * @see `readout::ROPDataContainer`
     *
     * This function operates as `makeROPdata() const`, except that copies
     * the specified value into all the entries of the container. Example:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * auto nHitsPerROP = geom->makeROPdata(0U);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template <typename T>
    readout::ROPDataContainer<T> makeROPdata(T const& defValue) const
    {
      return {Ncryostats(), MaxTPCsets(), MaxROPs(), defValue};
    }

    //
    // access
    //
    /// Returns whether we have the specified readout plane
    /// @return whether the readout plane is valid and exists
    bool HasROP(readout::ROPID const& ropid) const;

    /// Returns whether we have the specified readout plane
    /// @return whether the readout plane is valid and exists
    bool HasElement(readout::ROPID const& ropid) const { return HasROP(ropid); }

    //
    // mapping
    //
    /**
     * @brief Returns the ID of the ROP planeid belongs to
     * @param planeid ID of the wire plane
     * @return the ID of the ROP planeid belongs to
     *
     * If planeid is an invalid ID, an invalid ROP ID is returned.
     * If planeid is a valid ID (i.e. an ID whose isValid flag is set) that
     * points to a non-existent wire plane, the result is undefined.
     * Use HasPlaneID() to check if the wire plane actually exists.
     */
    readout::ROPID WirePlaneToROP(PlaneID const& planeid) const;

    /**
     * @brief Returns a list of ID of planes belonging to the specified ROP
     * @param ropid ID of the readout plane
     * @return list of ID of wire planes belonging to the specified ROP
     *
     * If ropid is an invalid ID, an empty list is returned.
     * If ropid is a valid ID (i.e. an ID whose isValid flag is set) that
     * points to a non-existent readout plane, the result is undefined.
     * Use HasROP() to check if the readout plane actually exists.
     */
    std::vector<PlaneID> ROPtoWirePlanes(readout::ROPID const& ropid) const;

    /**
     * @brief Returns a list of ID of TPCs the specified ROP spans
     * @param ropid ID of the readout plane
     * @return the list of TPC IDs, empty if readout plane ID is invalid
     *
     * Note that this check is performed on the validity of the readout plane
     * ID, that does not necessarily imply that the readout plane specified by
     * the ID actually exists. Check if the ROP exists with HasROP().
     * The behaviour on non-existing readout planes is undefined.
     */
    std::vector<TPCID> ROPtoTPCs(readout::ROPID const& ropid) const;

    /**
     * @brief Returns the ID of the first channel in the specified readout plane
     * @param ropid ID of the readout plane
     * @return ID of first channel, or raw::InvalidChannelID if ID is invalid
     *
     * Note that this check is performed on the validity of the readout plane
     * ID, that does not necessarily imply that the readout plane specified by
     * the ID actually exists. Check if the ROP exists with HasROP().
     * The behaviour for non-existing readout planes is undefined.
     */
    raw::ChannelID_t FirstChannelInROP(readout::ROPID const& ropid) const;

    ///
    /// iterators
    ///

    /// Initializes the specified ID with the ID of the first readout plane.
    void GetBeginID(readout::ROPID& id) const
    {
      GetBeginID(id.asTPCsetID());
      id.ROP = 0;
    }

    /// Initializes the specified ID with the invalid ID after the last ROP.
    void GetEndID(readout::ROPID& id) const
    {
      GetEndID(id.asTPCsetID());
      id.ROP = 0;
    }

    /// Sets the ID to the ID after the specified one.
    /// @return whether the ID is actually valid (validity flag is also set)
    bool IncrementID(readout::ROPID& id) const; // inline implementation

    /// Returns the ID of the first readout plane of the specified cryostat.
    readout::ROPID GetBeginROPID(CryostatID const& id) const { return {GetBeginTPCsetID(id), 0}; }

    /// Returns the (possibly invalid) ID after the last readout plane of the
    /// specified cryostat.
    readout::ROPID GetEndROPID(CryostatID const& id) const { return {GetEndTPCsetID(id), 0}; }

    /// Returns the ID of the first readout plane of the specified TPC set.
    readout::ROPID GetBeginROPID(readout::TPCsetID const& id) const { return {id, 0}; }

    /// Returns the (possibly invalid) ID after the last readout plane of the
    /// specified TPC set.
    readout::ROPID GetEndROPID(readout::TPCsetID const& id) const { return {GetNextID(id), 0}; }

    /**
     * @brief Returns the view of the channels in the specified readout plane
     * @param ropid readout plane ID
     * @return the type of signal on the specified ROP
     *
     * Returns the view (wire orientation) on the channels of specified readout
     * plane.
     * If ropid is an invalid ID, geo::kUnknown is returned.
     * If ropid is a valid ID (i.e. an ID whose isValid flag is set) that
     * points to a non-existent readout plane, the result is undefined.
     * Use HasROP() to check if the readout plane actually exists.
     */
    View_t View(readout::ROPID const& ropid) const;

    /**
     * @brief Returns the type of signal of channels in specified readout plane
     * @param ropid readout plane ID
     * @return the type of signal on the specified ROP
     *
     * Assumes that all the channels on the readout plane have the same signal
     * type.
     * If ropid is an invalid ID, geo::kMysteryType is returned.
     * If ropid is a valid ID (i.e. an ID whose isValid flag is set) that
     * points to a non-existent readout plane, the result is undefined.
     * Use HasROP() to check if the readout plane actually exists.
     */
    SigType_t SignalType(readout::ROPID const& ropid) const;

    /// @} Readout plane information

    /**
     * @name Optical readout channels
     * @anchor GeometryCoreOpDetChannel
     * @see @ref GeometryCoreOpDetGeometry "optical detector geometry information"
     */
    /// @todo add explanation of the different IDs
    /// @{

    //
    // group features
    //

    /// Number of electronics channels for all the optical detectors
    unsigned int NOpChannels() const;

    /// Largest optical channel number
    unsigned int MaxOpChannel() const;

    // Number of hardware channels for a given optical detector
    unsigned int NOpHardwareChannels(int opDet) const;

    //
    // access
    //

    /// Is this a valid OpChannel number?
    bool IsValidOpChannel(int opChannel) const;

    /// Convert detector number and hardware channel to unique channel
    unsigned int OpChannel(int detNum, int hardwareChannel) const;

    /// Convert unique channel to detector number
    unsigned int OpDetFromOpChannel(int opChannel) const;

    /// Convert unique channel to hardware channel
    unsigned int HardwareChannelFromOpChannel(int opChannel) const;

    /// Get unique opdet number from cryo and internal count
    unsigned int OpDetFromCryo(unsigned int o, unsigned int c) const;

    /// @} Optical readout channels

    //
    // unsorted methods
    //

    /**
     * @brief Returns whether a value is within the specified range
     * @param value the value to be tested
     * @param min the lower boundary
     * @param max the upper boundary
     * @return whether the value is within range
     *
     * If min is larger than max, they are swapped.
     * A tolerance of 10^-6 (absolute) is used.
     *
     * @todo Use wiggle instead of 10^-6
     * @todo resort source code for a bit of speed up
     */
    bool ValueInRange(double value, double min, double max) const
    {
      return unrelated::ValueInRange(value, min, max);
    }

    /// @name Geometry initialization
    /// @{

    /**
     * @brief Loads the geometry information from the specified files
     * @param gdmlfile path to file to be used for Geant4 simulation
     * @param rootfile path to file for internal geometry representation
     * @param builder algorithm to be used for the interpretation of geometry
     * @param bForceReload reload even if there is already a valid geometry
     * @see ApplyChannelMap()
     *
     * Both paths must directly resolve to an available file, as no search
     * is performed for them.
     *
     * The gdmlfile parameter does not have to necessarily be in GDML format,
     * as long as it's something supported by Geant4. This file is not used by
     * the geometry, but its path is provided on request by the simulation
     * modules (see LArSoft `LArG4` module).
     * The rootfile also does not need to be a ROOT file, but just anything
     * that TGeoManager::Import() supports. This file is parsed immediately
     * and the internal geometry representation is built out of it.
     *
     * @note After calling this method, the detector geometry information can
     * be considered complete, but the geometry service provider is not fully
     * initialized yet, since it's still necessary to provide or update the
     * channel mapping.
     */
    void LoadGeometryFile(std::string gdmlfile,
                          std::string rootfile,
                          GeometryBuilder& builder,
                          bool bForceReload = false);

    /**
     * @brief Loads the geometry information from the specified files
     * @param gdmlfile path to file to be used for Geant4 simulation
     * @param rootfile path to file for internal geometry representation
     * @param bForceReload reload even if there is already a valid geometry
     * @see ApplyChannelMap()
     *
     * This legacy version of `LoadGeometryFile()` uses a standard
     * `geo::GeometryBuilder` implementation.
     * Do not rely on it if you can avoid it.
     */
    void LoadGeometryFile(std::string gdmlfile, std::string rootfile, bool bForceReload = false);

    /**
     * @brief Initializes the geometry to work with this channel map
     * @param pChannelMap a pointer to the channel mapping algorithm to be used
     * @see LoadGeometryFile()
     *
     * The specified channel mapping is used with this geometry.
     * These modifications typically involve some resorting of the objects.
     *
     * This method needs to be called after LoadGeometryFile() to complete the
     * geometry initialization.
     */
    void ApplyChannelMap(std::unique_ptr<ChannelMapAlg> pChannelMap);
    /// @}

  private:
    /// Return the internal cryostat list
    CryostatList_t& Cryostats() { return fGeoData.cryostats; }
    CryostatList_t const& Cryostats() const { return fGeoData.cryostats; }

    /// Return the internal auxdet list
    AuxDetList_t& AuxDets() { return fGeoData.auxDets; }
    AuxDetList_t const& AuxDets() const { return fGeoData.auxDets; }

    GeometryData_t fGeoData; ///< The detector description data

    double fSurfaceY;          ///< The point where air meets earth for this detector.
    std::string fDetectorName; ///< Name of the detector.
    std::string fGDMLfile;     ///< path to geometry file used for Geant4 simulation
    std::string fROOTfile;     ///< path to geometry file for geometry in GeometryCore
    double fMinWireZDist;      ///< Minimum distance in Z from a point in which
                               ///< to look for the closest wire
    double fPositionWiggle;    ///< accounting for rounding errors when testing positions

    /// Configuration for the geometry builder
    /// (needed since builder is created after construction).
    fhicl::ParameterSet fBuilderParameters;
    std::unique_ptr<const ChannelMapAlg> fChannelMapAlg;
    ///< Object containing the channel to wire mapping

    // cached values
    std::set<View_t> allViews; ///< All views in the detector.

    std::vector<TGeoNode const*> FindDetectorEnclosure(
      std::string const& name = "volDetEnclosure") const;

    bool FindFirstVolume(std::string const& name, std::vector<const TGeoNode*>& path) const;

    /// Parses ROOT geometry nodes and builds LArSoft geometry representation.
    /// @param builder the algorithm to be used
    void BuildGeometry(GeometryBuilder& builder);

    /// Wire ID check for WireIDsIntersect methods
    bool WireIDIntersectionCheck(const WireID& wid1, const WireID& wid2) const;

    /// Runs the sorting of geometry with the sorter provided by channel mapping
    void SortGeometry(GeoObjectSorter const& sorter);

    /// Performs all the updates needed after sorting
    void UpdateAfterSorting();

    /// Deletes the detector geometry structures
    void ClearGeometry();

  }; // class GeometryCore

  /// @}
  // END Geometry group --------------------------------------------------------

} // namespace geo

template <typename Element, typename GEOIDITER>
geo::details::geometry_element_iterator<Element, GEOIDITER>::operator bool() const
{
  assert(geom);
  // FIXME: what if id_iter is invalid?
  return geom && geom->HasElement(*id_iter) && get() != nullptr;
}

template <typename Element, typename GEOIDITER>
auto geo::details::geometry_element_iterator<Element, GEOIDITER>::get() const -> ElementPtr_t
{
  assert(geom);
  // FIXME: what if id_iter is invalid?
  return geom->GetElementPtr(*id_iter);
}

//******************************************************************************
//*** inline implementation
//***
inline bool geo::GeometryCore::IncrementID(CryostatID& id) const
{
  ++id.Cryostat;
  if (id) id.isValid = HasCryostat(id); // if invalid already, it stays so
  return bool(id);
}

inline bool geo::GeometryCore::IncrementID(TPCID& id) const
{
  unsigned int const nTPCsInCryo = NTPC(id);
  if (++id.TPC < nTPCsInCryo) return bool(id); // if was invalid, it stays so
  // no more TPCs in this cryostat
  id.TPC = 0;
  return IncrementID(id.asCryostatID()); // also sets validity
}

inline bool geo::GeometryCore::IncrementID(PlaneID& id) const
{
  // this implementation is non-optimal, in that the cryostat lookup is
  // performed both here and, potentially, in IncrementID(TPCID)
  unsigned int const nPlanesInTPC = Nplanes(id);
  if (++id.Plane < nPlanesInTPC) return bool(id); // if was invalid, stays so
  // no more planes in this TPCs
  id.Plane = 0;
  return IncrementID(id.asTPCID()); // also sets validity
}

inline bool geo::GeometryCore::IncrementID(WireID& id) const
{
  // this implementation is non-optimal, in that the TPC lookup is
  // performed both here and, potentially, in IncrementID(PlaneID)
  unsigned int const nWiresInPlane = Nwires(id);
  if (++id.Wire < nWiresInPlane) return bool(id); // if was invalid, stays so
  // no more wires in this plane
  id.Wire = 0;
  return IncrementID(id.asPlaneID()); // also sets validity
}

inline bool geo::GeometryCore::IncrementID(readout::TPCsetID& id) const
{
  unsigned int const nTPCsetsInCryo = NTPCsets(id);
  if (++id.TPCset < nTPCsetsInCryo) return bool(id); // if was invalid, it stays so
  // no more TPC sets in this cryostat
  id.TPCset = 0;
  return IncrementID(id.asCryostatID()); // also sets validity
}

inline bool geo::GeometryCore::IncrementID(readout::ROPID& id) const
{
  // this implementation is non-optimal, in that the cryostat lookup is
  // performed both here and, potentially, in IncrementID(TPCsetID)
  unsigned int const nROPinTPC = NROPs(id);
  if (++id.ROP < nROPinTPC) return bool(id); // if was invalid, stays so
  // no more readout planes in this TPC set
  id.ROP = 0;
  return IncrementID(id.asTPCsetID()); // also sets validity
}

//******************************************************************************
//***  template implementation
//***
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// template member function specializations
namespace geo {

  // TPCID
  template <>
  inline TPCID GeometryCore::GetBeginID<TPCID, CryostatID>(CryostatID const& id) const
  {
    return GetBeginTPCID(id);
  }

  template <>
  inline TPCID GeometryCore::GetEndID<TPCID, CryostatID>(CryostatID const& id) const
  {
    return GetEndTPCID(id);
  }

  // PlaneID
  template <>
  inline PlaneID GeometryCore::GetBeginID<PlaneID, CryostatID>(CryostatID const& id) const
  {
    return GetBeginPlaneID(id);
  }

  template <>
  inline PlaneID GeometryCore::GetBeginID<PlaneID, TPCID>(TPCID const& id) const
  {
    return GetBeginPlaneID(id);
  }

  template <>
  inline PlaneID GeometryCore::GetEndID<PlaneID, CryostatID>(CryostatID const& id) const
  {
    return GetEndPlaneID(id);
  }

  template <>
  inline PlaneID GeometryCore::GetEndID<PlaneID, TPCID>(TPCID const& id) const
  {
    return GetEndPlaneID(id);
  }

  // WireID
  template <>
  inline WireID GeometryCore::GetBeginID<WireID, CryostatID>(CryostatID const& id) const
  {
    return GetBeginWireID(id);
  }

  template <>
  inline WireID GeometryCore::GetBeginID<WireID, TPCID>(TPCID const& id) const
  {
    return GetBeginWireID(id);
  }

  template <>
  inline WireID GeometryCore::GetBeginID<WireID, PlaneID>(PlaneID const& id) const
  {
    return GetBeginWireID(id);
  }

  template <>
  inline WireID GeometryCore::GetEndID<WireID, CryostatID>(CryostatID const& id) const
  {
    return GetEndWireID(id);
  }

  template <>
  inline WireID GeometryCore::GetEndID<WireID, TPCID>(TPCID const& id) const
  {
    return GetEndWireID(id);
  }

  template <>
  inline WireID GeometryCore::GetEndID<WireID, PlaneID>(PlaneID const& id) const
  {
    return GetEndWireID(id);
  }

  // TPCsetID
  template <>
  inline readout::TPCsetID GeometryCore::GetBeginID<readout::TPCsetID, CryostatID>(
    CryostatID const& id) const
  {
    return GetBeginTPCsetID(id);
  }

  template <>
  inline readout::TPCsetID GeometryCore::GetEndID<readout::TPCsetID, CryostatID>(
    CryostatID const& id) const
  {
    return GetEndTPCsetID(id);
  }

  // ROPID
  template <>
  inline readout::ROPID GeometryCore::GetBeginID<readout::ROPID, CryostatID>(
    CryostatID const& id) const
  {
    return GetBeginROPID(id);
  }

  template <>
  inline readout::ROPID GeometryCore::GetEndID<readout::ROPID, CryostatID>(
    CryostatID const& id) const
  {
    return GetEndROPID(id);
  }

  template <>
  inline readout::ROPID GeometryCore::GetBeginID<readout::ROPID, readout::TPCsetID>(
    readout::TPCsetID const& id) const
  {
    return GetBeginROPID(id);
  }

  template <>
  inline readout::ROPID GeometryCore::GetEndID<readout::ROPID, readout::TPCsetID>(
    readout::TPCsetID const& id) const
  {
    return GetEndROPID(id);
  }

} // namespace geo

template <typename Point>
geo::GeometryCore::Segment<Point> geo::GeometryCore::WireEndPoints(WireID const& wireid) const
{
  WireGeo const& wire = Wire(wireid);
  return {wire.GetStart<Point>(), wire.GetEnd<Point>()};
}

//------------------------------------------------------------------------------
template <typename Stream>
void geo::GeometryCore::Print(Stream&& out, std::string indent /* = "  " */) const
{

  out << "Detector " << DetectorName() << " has " << Ncryostats() << " cryostats and " << NAuxDets()
      << " auxiliary detectors:";

  auto const& detEnclosureBox = DetectorEnclosureBox();
  out << "\n"
      << indent << "Detector enclosure: " << detEnclosureBox.Min() << " -- "
      << detEnclosureBox.Max() << " cm => ( " << detEnclosureBox.SizeX() << " x "
      << detEnclosureBox.SizeY() << " x " << detEnclosureBox.SizeZ() << " ) cm^3";

  for (auto const& cryostat : Iterate<CryostatGeo>()) {
    out << "\n" << indent;
    cryostat.PrintCryostatInfo(std::forward<Stream>(out), indent + "  ", cryostat.MaxVerbosity);

    const unsigned int nTPCs = cryostat.NTPC();
    for (unsigned int t = 0; t < nTPCs; ++t) {
      const TPCGeo& tpc = cryostat.TPC(t);

      out << "\n" << indent << "  ";
      tpc.PrintTPCInfo(std::forward<Stream>(out), indent + "    ", tpc.MaxVerbosity);

      const unsigned int nPlanes = tpc.Nplanes();
      for (unsigned int p = 0; p < nPlanes; ++p) {
        const PlaneGeo& plane = tpc.Plane(p);
        const unsigned int nWires = plane.Nwires();

        out << "\n" << indent << "    ";
        plane.PrintPlaneInfo(std::forward<Stream>(out), indent + "      ", plane.MaxVerbosity);
        SigType_t const sigType = SignalType(plane.ID());
        out << "\n"
            << indent << "      "
            << "signal type: " << SignalTypeName(sigType) << " (" << static_cast<int>(sigType)
            << ")";

        for (unsigned int w = 0; w < nWires; ++w) {
          const WireGeo& wire = plane.Wire(w);
          WireID wireID(plane.ID(), w);

          // the wire should be aligned on z axis, half on each side of 0,
          // in its local frame
          out << "\n" << indent << "      " << wireID << " ";
          wire.PrintWireInfo(std::forward<Stream>(out), indent + "      ", wire.MaxVerbosity);
        } // for wire
      }   // for plane
    }     // for TPC

    unsigned int nOpDets = cryostat.NOpDet();
    for (unsigned int iOpDet = 0; iOpDet < nOpDets; ++iOpDet) {
      OpDetGeo const& opDet = cryostat.OpDet(iOpDet);
      out << "\n" << indent << "  [OpDet #" << iOpDet << "] ";
      opDet.PrintOpDetInfo(std::forward<Stream>(out), indent + "  ", opDet.MaxVerbosity);
    } // for
  }   // for cryostat

  unsigned int const nAuxDets = NAuxDets();
  for (unsigned int iDet = 0; iDet < nAuxDets; ++iDet) {
    AuxDetGeo const& auxDet = AuxDet(iDet);

    out << "\n" << indent << "[#" << iDet << "] ";
    auxDet.PrintAuxDetInfo(std::forward<Stream>(out), indent + "  ", auxDet.MaxVerbosity);

    unsigned int const nSensitive = auxDet.NSensitiveVolume();
    switch (nSensitive) {
    case 0: break;
    case 1: {
      AuxDetSensitiveGeo const& auxDetS = auxDet.SensitiveVolume(0U);
      out << "\n" << indent << "  ";
      auxDetS.PrintAuxDetInfo(std::forward<Stream>(out), indent + "    ", auxDetS.MaxVerbosity);
      break;
    }
    default:
      for (unsigned int iSens = 0; iSens < nSensitive; ++iSens) {
        out << "\n" << indent << "[#" << iSens << "] ";
        AuxDetSensitiveGeo const& auxDetS = auxDet.SensitiveVolume(iSens);
        auxDetS.PrintAuxDetInfo(std::forward<Stream>(out), indent + "    ", auxDetS.MaxVerbosity);
      } // for
      break;
    } // if sensitive detectors

  } // for auxiliary detector

  out << '\n';

} // geo::GeometryCore::Print()

#endif // LARCOREALG_GEOMETRY_GEOMETRYCORE_H
