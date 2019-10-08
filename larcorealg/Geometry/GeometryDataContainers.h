/**
 * @file   larcorealg/Geometry/GeometryDataContainers.h
 * @brief  Containers to hold one datum per TPC or plane.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   January 2nd, 2018
 * @ingroup Geometry
 *
 * This is a header-only library.
 */

#ifndef LARCOREALG_GEOMETRY_GEOMETRYDATACONTAINERS_H
#define LARCOREALG_GEOMETRY_GEOMETRYDATACONTAINERS_H

// LArSoft libraries
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"

// C/C++ standard libraries
#include <vector>
#include <array>
#include <initializer_list>
#include <string>
#include <utility> // std::forward()
#include <algorithm> // std::fill(), std::for_each()
#include <stdexcept> // std::out_of_range
#include <cassert>


namespace geo {
  
  template <typename T, typename IDType>
  class GeoIDdataContainer;
  
  template <typename T>
  class TPCDataContainer;
  
  template <typename T>
  class PlaneDataContainer;
  
  // ---------------------------------------------------------------------------
  namespace details {
    
    template <typename T>
    class GeoContainerData;
    
    /// Returns a STL array of size `N` filled with `values` from the argument.
    template <std::size_t N, typename T>
    auto initializerListToArray(std::initializer_list<T> values);
    
  } // namespace details
  // ---------------------------------------------------------------------------
  
} // namespace geo


// --- BEGIN Geometry data containers ----------------------------------------

/// @name Geometry data containers
/// @ingroup Geometry
/// @{

/** **************************************************************************
 * @brief Container with one element per geometry TPC.
 * @tparam T type of the contained datum
 * @see `geo::GeometryCore::makeTPCData`
 *
 * The container is of fixed size and can't be neither resized nor freed
 * before destruction.
 *
 * This example creates a "map" of tracks starting on each TPC:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * auto const* geom = lar::providerFrom<geo::GeometryCore>();
 * geo::TPCDataContainer<std::vector<recob::Track const*>> TracksPerTPC
 *   (geom->NCryostats(), geom->MaxTPCs());
 *
 * for (recob::Track const& track: tracks) {
 *   geo::TPCGeo const* tpc = geom->PositionToTPCptr(track.Start());
 *   if (tpc) TracksPerTPC[tpc->ID()].push_back(tpc);
 * } // for
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *
 * Assumptions
 * ============
 *
 * The following assumptions should be considered unchecked, and the behavior
 * when they are violated undefined (but note that in debug mode some of them
 * might be actually checked):
 * * the container assumes the same number of TPCs in each cryostat. While
 *   this is not effectively a necessary condition, keep in mind that this
 *   container has no notion whether a given TPC actually exists in the
 *   geometry or not
 * * at least one element is expected to be present
 *
 */
template <typename T, typename IDType>
class geo::GeoIDdataContainer {
  
  /// Type of data container helper.
  using Container_t = details::GeoContainerData<T>;
  
  using ID_t = IDType; ///< Type used as ID for this container.
  
    public:
  
  /// @{
  /// @name STL container types.

  // TODO we don't offer iterators so far

  using value_type             = typename Container_t::value_type            ;
  using reference              = typename Container_t::reference             ;
  using const_reference        = typename Container_t::const_reference       ;
  using pointer                = typename Container_t::pointer               ;
  using const_pointer          = typename Container_t::const_pointer         ;
//     using iterator               = typename Container_t::iterator              ;
//     using const_iterator         = typename Container_t::const_iterator        ;
//     using reverse_iterator       = typename Container_t::reverse_iterator      ;
//     using const_reverse_iterator = typename Container_t::const_reverse_iterator;
  using difference_type        = typename Container_t::difference_type       ;
  using size_type              = typename Container_t::size_type             ;

  /// @}

  /**
   * @brief Prepares the container with default-constructed data.
   * @param dims number of elements on all levels of the container
   * 
   * The size of each dimension is specified by the corresponding number,
   * starting from the size of the outer dimension (cryostat).
   *
   * The container is sized to host data for all the elements.
   * Each element in the container is default-constructed.
   */
  GeoIDdataContainer(std::initializer_list<unsigned int> dims);

  /**
   * @brief Prepares the container with default-constructed data.
   * @param dims number of elements on all levels of the container
   * 
   * The size of each dimension is specified by the corresponding number,
   * starting from the size of the outer dimension (cryostat).
   *
   * The container is sized to host data for all the elements.
   * Each element in the container is default-constructed.
   */
  GeoIDdataContainer(
    std::initializer_list<unsigned int> dims,
    value_type const& defValue
    );


  // --- BEGIN Container status query ----------------------------------------
  /// @name Container status query
  /// @{

  /// Returns the number of elements in the container.
  size_type size() const;

  /// Returns the number of elements the container has memory for.
  size_type capacity() const;

  /// Returns whether the container has no elements (`false` by assumptions).
  bool empty() const;

  /// Dimensions of the `Level` dimension of this container.
  template <std::size_t Level>
  unsigned int dimSize();
  
  /// Dimensions of the ID of this container.
  static constexpr unsigned int dimensions();
  
  /// Returns whether this container hosts data for the specified ID.
  template <typename GeoID>
  bool hasElement(GeoID const& id) const;

  /// Returns the ID of the first element with GeoID type.
  template <typename GeoID = ID_t>
  GeoID firstID() const;

  /// Returns the ID of the last covered element with GeoID type.
  template <typename GeoID = ID_t>
  GeoID lastID() const;

  /// @}
  // --- END Container status query ------------------------------------------


  // --- BEGIN Element access ------------------------------------------------
  /// @name Element access
  /// @{

  /// Returns the element for the specified geometry element.
  
  reference operator[](ID_t const& id);

  /// Returns the element for the specified geometry element (read-only).
  const_reference operator[](ID_t const& id) const;

  /// Returns the element for the specified geometry element.
  /// @throw std::out_of_range if element `id` is not within the container range
  reference at(ID_t const& id);

  /// Returns the element for the specified geometry element (read-only).
  /// @throw std::out_of_range if element `id` is not within the container range
  const_reference at (ID_t const& id) const;


  /// Returns the element for the first ID (unchecked).
  reference first();

  /// Returns the element for the first ID (unchecked).
  const_reference first() const;


  /// Returns the element for the last ID (unchecked).
  reference last();

  /// Returns the element for the last ID (unchecked).
  const_reference last() const;

  /// @}
  // --- END Element access --------------------------------------------------


  // --- BEGIN Data modification ---------------------------------------------
  /**
   * @name Data modification
   * 
   * In general, each single element can be accessed and changed.
   * In addition, this section includes methods acting on multiple elements
   * at once.
   */
  /// @{
  
  /// Sets all elements to the specified `value` (copied).
  void fill(value_type value);
  
  /// Sets all the elements to a default-constructed `value_type`.
  void clear();
  
  
  /**
   * @brief Applies an operation on all elements.
   * @tparam Op type of operation
   * @param op Operation
   * @return the operation object after operations took place
   * 
   * The operation `op` is a unary functor, i.e. an object that supports the
   * call to `op(value_type&)`.
   * 
   * The return values of `op` calls are discarded.
   */
  template <typename Op>
  Op apply(Op&& op);
  
  /**
   * @brief Applies an operation on all elements.
   * @tparam Op type of operation
   * @param op Operation
   * @return the operation object after operations took place
   * 
   * The operation `op` is a unary functor, i.e. an object that supports the
   * call to `op(value_type const&)`.
   * 
   * The return values of `op` calls are discarded.
   * Note that while the elements of this container can't be modified, the
   * operation itself still can if not constant, and it is returned.
   */
  template <typename Op>
  decltype(auto) apply(Op&& op) const;
  
  /// @}
  // --- END Element access --------------------------------------------------


    private:
  ///< Type of dimension sizes.
  using Dimensions_t = std::array<unsigned int, dimensions()>;
  
  Dimensions_t fN; ///< Number of maximum entries per ID level.
  
  Container_t fData; ///< Data storage.


  /// Returns the internal index of the specified TPC in the storage area.
  size_type index(ID_t const& id) const;
  
  template <std::size_t Level, typename GeoID>
  size_type indexLevel(GeoID const& id) const;
  
  /// Returns whether all levels of `id` up to `Level` are within range.
  template <std::size_t Level, typename GeoID>
  bool hasElementLevel(GeoID const& id) const;

  /// Computes the expected size of this container.
  size_type computeSize() const;

  /// Returns the number of elements at the specified `Level`.
  /// @param dimSizes the sizes of each of the levels
  template <std::size_t Level, typename Dims>
  static size_type sizeLevel(Dims const& dimSizes);
  
}; // class geo::GeoIDdataContainer<>


//------------------------------------------------------------------------------
/** 
 * @brief Container with one element per geometry TPC.
 * @tparam T type of the contained datum
 * @see `geo::GeometryCore::makeTPCData`
 *
 * The container is of fixed size and can't be neither resized nor freed
 * before destruction.
 *
 * This example creates a "map" of tracks starting on each TPC:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * auto const* geom = lar::providerFrom<geo::GeometryCore>();
 * geo::TPCDataContainer<std::vector<recob::Track const*>> TracksPerTPC
 *   (geom->NCryostats(), geom->MaxTPCs());
 *
 * for (recob::Track const& track: tracks) {
 *   geo::TPCGeo const* tpc = geom->PositionToTPCptr(track.Start());
 *   if (tpc) TracksPerTPC[tpc->ID()].push_back(tpc);
 * } // for
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *
 * Assumptions
 * ============
 *
 * The following assumptions should be considered unchecked, and the behavior
 * when they are violated undefined (but note that in debug mode some of them
 * might be actually checked):
 * * the container assumes the same number of TPCs in each cryostat. While
 *   this is not effectively a necessary condition, keep in mind that this
 *   container has no notion whether a given TPC actually exists in the
 *   geometry or not
 * * at least one element is expected to be present
 *
 */
template <typename T>
class geo::TPCDataContainer: public geo::GeoIDdataContainer<T, geo::TPCID> {
  
  using BaseContainer_t = geo::GeoIDdataContainer<T, geo::TPCID>;
  
    public:
  
  using value_type = typename BaseContainer_t::value_type;
  
  /**
   * @brief Prepares the container with default-constructed data.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs
   *
   * The container is sized to host data for `nCryo` cryostats, each with
   * `nTPCs` TPCs. Each element in the container is default-constructed.
   */
  TPCDataContainer(unsigned int nCryo, unsigned int nTPCs)
    : BaseContainer_t({ nCryo, nTPCs })
    {}

  /**
   * @brief Prepares the container with copies of the specified default value.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs
   * @param defValue the value to be replicated
   *
   * The container is sized to host data for `nCryo` cryostats, each with
   * `nTPCs` TPCs. Each element in the container is a copy of defValue.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const* geom = lar::providerFrom<geo::GeometryCore>();
   * geo::TPCDataContainer<unsigned int> PlanesPerTPC
   *   (geom->NCryostats(), geom->MaxTPCs(), 3U);
   * for (geo::TPCGeo const& TPC: geom->IterateTPC())
   *   assert(PlanesPerTPC[TPC.ID()] == TPC.Nplanes());
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  TPCDataContainer
    (unsigned int nCryo, unsigned int nTPCs, value_type const& defValue)
    : BaseContainer_t({ nCryo, nTPCs }, defValue)
    {}


  // --- BEGIN Container status query ------------------------------------------
  /// @name Container status query
  /// @{

  /// Returns whether this container hosts data for the specified cryostat.
  bool hasCryostat(geo::CryostatID const& cryoid) const
    { return BaseContainer_t::hasElement(cryoid); }

  /// Returns whether this container hosts data for the specified TPC.
  bool hasTPC(geo::TPCID const& tpcid) const
    { return BaseContainer_t::hasElement(tpcid); }

  /// @}
  // --- END Container status query --------------------------------------------


}; // class geo::details::TPCDataContainer<>


//------------------------------------------------------------------------------
/** 
 * @brief Container with one element per geometry wire plane.
 * @tparam T type of the contained datum
 * @see `geo::GeometryCore::makePlaneData`
 *
 * The container is of fixed size and can't be neither resized nor freed
 * before destruction.
 *
 * This example creates a "map" of hits on each wire plane:
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * auto const* geom = lar::providerFrom<geo::GeometryCore>();
 * geo::PlaneDataContainer<std::vector<recob::Hit const*>> hitsPerPlane
 *   (geom->NCryostats(), geom->MaxTPCs(), geom->MaxPlanes());
 *
 * for (recob::Hit const& hit: hits) {
 *   if (hit.WireID()) hitsPerPlane[hit.WireID().planeID()].push_back(&hit);
 * } // for
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 *
 * Assumptions
 * ============
 *
 * The following assumptions should be considered unchecked, and the behavior
 * when they are violated undefined (but note that in debug mode some of them
 * might be actually checked):
 * * the container assumes the same number of planes in all TPCs, and of TPCs
 *   on all cryostats. While this is not effectively a necessary condition,
 *   keep in mind that this container has no notion whether a given plane or
 *   TPC actually exists in the geometry or not
 * * at least one element is expected to be present
 *
 */
template <typename T>
class geo::PlaneDataContainer: public geo::GeoIDdataContainer<T, geo::PlaneID>
{

  /// Base class.
  using BaseContainer_t = geo::GeoIDdataContainer<T, geo::PlaneID>;
  
    public:

  /**
   * @brief Prepares the container with default-constructed data.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs per cryostat
   * @param nPlanes number of planes per TPC
   *
   * The container is sized to host data for `nCryo` cryostats, each with
   * `nTPCs` TPCs, each one with `nPlanes` wire planes. Each element in the
   * container is default-constructed.
   */
  PlaneDataContainer
    (unsigned int nCryo, unsigned int nTPCs, unsigned int nPlanes)
    : BaseContainer_t({ nCryo, nTPCs, nPlanes })
    {}

  /**
   * @brief Prepares the container with copies of the specified default value.
   * @param nCryo number of cryostats
   * @param nTPCs number of TPCs
   * @param nPlanes number of planes per TPC
   * @param defValue the value to be replicated
   *
   * The container is sized to host data for `nCryo` cryostats, each with
   * `nTPCs` TPCs, and each of them with `nPlanes` planes. Each element in
   * the container is a copy of `defValue`.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const* geom = lar::providerFrom<geo::GeometryCore>();
   * geo::TPCDataContainer<unsigned int> countPerPlane
   *   (geom->NCryostats(), geom->MaxTPCs(), geom->MaxPlanes(), 0U);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  PlaneDataContainer(
    unsigned int nCryo, unsigned int nTPCs, unsigned int nPlanes,
    T const& defValue
    )
    : BaseContainer_t{ { nCryo, nTPCs, nPlanes }, defValue }
    {}


  // --- BEGIN Container status query ------------------------------------------
  /// @name Container status query
  /// @{

  /// Returns whether this container hosts data for the specified cryostat.
  bool hasCryostat(geo::CryostatID const& cryoid) const
    { return BaseContainer_t::hasElement(cryoid); }

  /// Returns whether this container hosts data for the specified TPC.
  bool hasTPC(geo::TPCID const& tpcid) const
    { return BaseContainer_t::hasElement(tpcid); }

  /// Returns whether this container hosts data for the specified plane.
  bool hasPlane(geo::PlaneID const& planeid) const
    { return BaseContainer_t::hasElement(planeid); }

  /// @}
  // --- END Container status query --------------------------------------------


}; // class geo::PlaneDataContainer


/// @}
// --- END Geometry data containers --------------------------------------------
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
//--- Template implementation
//------------------------------------------------------------------------------
template <std::size_t N, typename T>
auto geo::details::initializerListToArray(std::initializer_list<T> values) {
  std::array<T, N> data;
  std::copy(values.begin(), values.end(), data.begin());
  return data;
} // geo::details::initializerListToArray()


//------------------------------------------------------------------------------
//--- geo::details::GeoContainerData
//------------------------------------------------------------------------------
template <typename T>
class geo::details::GeoContainerData {

  using Container_t = std::vector<T>;

    public:

  // --- BEGIN STL container types ---------------------------------------------
  /// @{
  /// @name STL container types.

  using value_type             = typename Container_t::value_type            ;
  using reference              = typename Container_t::reference             ;
  using const_reference        = typename Container_t::const_reference       ;
  using pointer                = typename Container_t::pointer               ;
  using const_pointer          = typename Container_t::const_pointer         ;
  using iterator               = typename Container_t::iterator              ;
  using const_iterator         = typename Container_t::const_iterator        ;
  using reverse_iterator       = typename Container_t::reverse_iterator      ;
  using const_reverse_iterator = typename Container_t::const_reverse_iterator;
  using difference_type        = typename Container_t::difference_type       ;
  using size_type              = typename Container_t::size_type             ;

  /// @}
  // --- END STL container types -----------------------------------------------
  

  using index_type = size_type; ///< Type used internally (so far) for indexing.
  
  // --- BEGIN Constructors ----------------------------------------------------
  
  /// Prepares the container with default-constructed data.
  GeoContainerData(size_type size): fData(size) {}

  // Prepares the container with copies of the specified default value.
  GeoContainerData(size_type size, value_type const& defValue)
    : fData(size, defValue) {}
  
  // --- END Constructors ------------------------------------------------------
  
  
  // --- BEGIN Container status query ------------------------------------------
  /// @name Container status query
  /// @{

  /// Returns the number of elements in the container.
  size_type size() const { return fData.size(); }

  /// Returns the number of elements the container has memory for.
  size_type capacity() const { return fData.capacity(); }

  /// Returns whether the container has no elements (`false` by assumptions).
  bool empty() const { return fData.empty(); }

  /// @}
  // --- END Container status query --------------------------------------------
  
  
  // --- BEGIN Data modification -----------------------------------------------
  /**
   * @name Data modification
   * 
   * In general, each single element can be accessed and changed.
   * In addition, this section includes methods acting on multiple elements
   * at once.
   */
  /// @{
  
  /// Sets all elements to the specified `value` (copied).
  void fill(value_type value)
    { std::fill(fData.begin(), fData.end(), value); }
  
  /// Sets all the elements to a default-constructed `value_type`.
  void clear() { fill(value_type{}); }
  
  /**
   * @brief Applies an operation on all elements.
   * @tparam Op type of operation
   * @param op Operation
   * @return the operation object after operations took place
   * 
   * The operation `op` is a unary functor, i.e. an object that supports the
   * call to `op(value_type&)`.
   * 
   * The return values of `op` calls are discarded.
   */
  template <typename Op>
  Op apply(Op&& op)
    { for (auto& data: fData) op(data); return op; }
  
  /**
   * @brief Applies an operation on all elements.
   * @tparam Op type of operation
   * @param op Operation
   * @return the operation object after operations took place
   * 
   * The operation `op` is a unary functor, i.e. an object that supports the
   * call to `op(value_type const&)`.
   * 
   * The return values of `op` calls are discarded.
   */
  template <typename Op>
  Op apply(Op&& op) const
    { for (auto const& data: fData) op(data); return op; }
  
  
  /// @}
  // --- END Element access ----------------------------------------------------


  // --- BEGIN Element access --------------------------------------------------
  /// @name Element access
  /// @{
  
  /// Returns the element for the specified index.
  reference operator[](index_type index) { return fData[index]; }
  
  /// Returns the element for the specified index (read-only).
  const_reference operator[](index_type index) const { return fData[index]; }
  
  /// @}
  // --- END Element access ----------------------------------------------------


  /// Returns whether the specified value is between `0` and the upper limit.
  template <typename Value, typename Upper>
  static bool bounded(Value v, Upper upper)
    { return (v >= 0) && (static_cast<size_type>(v) < upper); }
  
    private:

  Container_t fData; ///< Data storage area.

}; // class geo::details::GeoContainerData


//------------------------------------------------------------------------------
//--- geo::GeoIDdataContainer
//------------------------------------------------------------------------------
template <typename T, typename IDType>
geo::GeoIDdataContainer<T, IDType>::GeoIDdataContainer
  (std::initializer_list<unsigned int> dims)
  : fN(details::initializerListToArray<dimensions()>(dims))
  , fData(computeSize())
{
  assert(dims.size() == dimensions()); // can't be static
  assert(!fData.empty()); 
}


//------------------------------------------------------------------------------
template <typename T, typename IDType>
geo::GeoIDdataContainer<T, IDType>::GeoIDdataContainer(
  std::initializer_list<unsigned int> dims,
  value_type const& defValue
  )
  : fN(details::initializerListToArray<dimensions()>(dims))
  , fData(computeSize(), defValue)
{
  assert(dims.size() == dimensions()); // can't be static
  assert(!fData.empty()); 
}


//------------------------------------------------------------------------------
template <typename T, typename IDType>
auto geo::GeoIDdataContainer<T, IDType>::size() const -> size_type 
  { return fData.size(); }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
auto geo::GeoIDdataContainer<T, IDType>::capacity() const -> size_type 
  { return fData.capacity(); }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
bool geo::GeoIDdataContainer<T, IDType>::empty() const
  { return fData.empty(); }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
template <std::size_t Level>
unsigned int geo::GeoIDdataContainer<T, IDType>::dimSize() {
  if constexpr (Level >= dimensions()) return 0U; // technically it would be 1...
  else return fN[Level];
} // geo::GeoIDdataContainer<>::dimSize()


//------------------------------------------------------------------------------
template <typename T, typename IDType>
constexpr unsigned int geo::GeoIDdataContainer<T, IDType>::dimensions()
  { return IDType::Level + 1; }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
template <typename GeoID>
bool geo::GeoIDdataContainer<T, IDType>::hasElement(GeoID const& id) const
  { return hasElementLevel<GeoID::Level>(id); }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
template <typename GeoID /* = ID_t */>
GeoID geo::GeoIDdataContainer<T, IDType>::firstID() const {
  if constexpr (GeoID::Level == 0) return GeoID(0U);
  else return GeoID(firstID<typename GeoID::ParentID_t>(), 0U);
} // geo::GeoIDdataContainer<>::firstID()


//------------------------------------------------------------------------------
template <typename T, typename IDType>
template <typename GeoID /* = ID_t */>
GeoID geo::GeoIDdataContainer<T, IDType>::lastID() const {
  if constexpr (GeoID::Level == 0) return GeoID(fN[GeoID::Level] - 1U);
  else
    return GeoID(lastID<typename GeoID::ParentID_t>(), fN[GeoID::Level] - 1U);
} // geo::GeoIDdataContainer<>::lastID()


//------------------------------------------------------------------------------
template <typename T, typename IDType>
auto geo::GeoIDdataContainer<T, IDType>::operator[](ID_t const& id) -> reference
  { return fData[index(id)]; }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
auto geo::GeoIDdataContainer<T, IDType>::operator[](ID_t const& id) const
  -> const_reference
  { return fData[index(id)]; }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
auto geo::GeoIDdataContainer<T, IDType>::at(ID_t const& id) -> reference {
  if (hasElement(id)) return operator[](id);
  throw std::out_of_range("No data for " + std::string(id));
} // geo::GeoIDdataContainer<>::at()


//------------------------------------------------------------------------------
template <typename T, typename IDType>
auto geo::GeoIDdataContainer<T, IDType>::at(ID_t const& id) const
  -> const_reference
{
  if (hasElement(id)) return operator[](id);
  throw std::out_of_range("No data for " + std::string(id));
} // geo::GeoIDdataContainer<>::at() const


//------------------------------------------------------------------------------
template <typename T, typename IDType>
auto geo::GeoIDdataContainer<T, IDType>::first() -> reference
  { return operator[](firstID()); }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
auto geo::GeoIDdataContainer<T, IDType>::first() const -> const_reference
  { return operator[](firstID()); }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
auto geo::GeoIDdataContainer<T, IDType>::last() -> reference
  { return operator[](lastID()); }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
auto geo::GeoIDdataContainer<T, IDType>::last() const -> const_reference
  { return operator[](lastID()); }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
void geo::GeoIDdataContainer<T, IDType>::fill(value_type value)
  { fData.fill(value); }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
void geo::GeoIDdataContainer<T, IDType>::clear()
  { fData.clear(); }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
template <typename Op>
Op geo::GeoIDdataContainer<T, IDType>::apply(Op&& op)
  { return fData.apply(std::forward<Op>(op)); }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
template <typename Op>
decltype(auto) geo::GeoIDdataContainer<T, IDType>::apply(Op&& op) const
  { return fData.apply(std::forward<Op>(op)); }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
auto geo::GeoIDdataContainer<T, IDType>::index(ID_t const& id) const
  -> size_type
  { return indexLevel<ID_t::Level>(id); }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
template <std::size_t Level, typename GeoID>
auto geo::GeoIDdataContainer<T, IDType>::indexLevel(GeoID const& id) const
  -> size_type
{
  if constexpr (Level == 0) return id.template getIndex<0U>();
  else {
    return
      indexLevel<(Level-1U)>(id) * fN[Level] + id.template getIndex<Level>();
  }
} // geo::GeoIDdataContainer<>::indexLevel()


//------------------------------------------------------------------------------
template <typename T, typename IDType>
template <std::size_t Level, typename GeoID>
bool geo::GeoIDdataContainer<T, IDType>::hasElementLevel(GeoID const& id) const
{
  if (!Container_t::bounded(id.template getIndex<Level>(), fN[Level]))
    return false;
  if constexpr (Level == 0U) return true;
  else return hasElementLevel<(Level-1U)>(id);
} // geo::GeoIDdataContainer<>::hasElementLevel()


//------------------------------------------------------------------------------
template <typename T, typename IDType>
auto geo::GeoIDdataContainer<T, IDType>::computeSize() const -> size_type
  { return sizeLevel<0U>(fN); }


//------------------------------------------------------------------------------
template <typename T, typename IDType>
template <std::size_t Level, typename Dims>
auto geo::GeoIDdataContainer<T, IDType>::sizeLevel(Dims const& dimSizes)
  -> size_type
{
  if constexpr (Level >= dimensions()) return 1U;
  else return sizeLevel<(Level+1U)>(dimSizes) * dimSizes[Level];
} // geo::GeoIDdataContainer<>::sizeLevel()


//------------------------------------------------------------------------------


#endif // LARCOREALG_GEOMETRY_GEOMETRYDATACONTAINERS_H
