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
#include <string>
#include <utility> // std::forward()
#include <algorithm> // std::fill(), std::for_each()
#include <stdexcept> // std::out_of_range
#include <cassert>


namespace geo {
  
  namespace details {
    
    template <typename T>
    class GeoContainerData;
    
  } // namespace details
  
  
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
  template <typename T>
  class TPCDataContainer {
    
    /// Type of data container helper.
    using Container_t = details::GeoContainerData<T>;
    
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
     * @param nCryo number of cryostats
     * @param nTPCs number of TPCs
     *
     * The container is sized to host data for `nCryo` cryostats, each with
     * `nTPCs` TPCs. Each element in the container is default-constructed.
     */
    TPCDataContainer(unsigned int nCryo, unsigned int nTPCs)
      : fNCryo(nCryo)
      , fNTPCs(nTPCs)
      , fData(computeSize())
      { assert(!fData.empty()); }

    /**
     * @brief Prepares the container with copies of the specified default value.
     * @param nCryo number of cryostats
     * @param nTPCs number of TPCs
     * @param defValue the value to be replicated
     *
     * The container is sized to host data for `nCryo` cryostats, each with
     * `nTPCs` TPCs. Each element in the container is a copy of defValue.
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * geo::TPCDataContainer<unsigned int> PlanesPerTPC
     *   (geom->NCryostats(), geom->MaxTPCs(), 3U);
     * for (geo::TPCGeo const& TPC: geom->IterateTPC())
     *   assert(PlanesPerTPC[TPC.ID()] == TPC.Nplanes());
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    TPCDataContainer
      (unsigned int nCryo, unsigned int nTPCs, value_type const& defValue)
      : fNCryo(nCryo)
      , fNTPCs(nTPCs)
      , fData(computeSize(), defValue)
      { assert(!fData.empty()); }


    // --- BEGIN Container status query ----------------------------------------
    /// @name Container status query
    /// @{

    /// Returns the number of elements in the container.
    size_type size() const { return fData.size(); }

    /// Returns the number of elements the container has memory for.
    size_type capacity() const { return fData.capacity(); }

    /// Returns whether the container has no elements (`false` by assumptions).
    bool empty() const { return fData.empty(); }

    /// Returns whether this container hosts data for the specified cryostat.
    bool hasCryostat(geo::CryostatID const& cryoid) const
      { return Container_t::bounded(cryoid.Cryostat, fNCryo); }

    /// Returns whether this container hosts data for the specified TPC.
    bool hasTPC(geo::TPCID const& id) const
      { return hasCryostat(id) && Container_t::bounded(id.TPC, fNTPCs); }

    /// Returns the ID of the first TPC.
    geo::TPCID firstID() const { return { 0U, 0U }; }

    /// Returns the ID of the last covered TPC.
    geo::TPCID lastID() const { return { fNCryo - 1U, fNTPCs - 1U }; }

    /// @}
    // --- END Container status query ------------------------------------------


    // --- BEGIN Element access ------------------------------------------------
    /// @name Element access
    /// @{

    /// Returns the element for the specified TPC.
    reference operator[](geo::TPCID const& id)
      { return fData[index(id)]; }

    /// Returns the element for the specified TPC (read-only).
    const_reference operator[](geo::TPCID const& id) const
      { return fData[index(id)]; }

    /// Returns the element for the specified TPC.
    /// @throw std::out_of_range if the TPC id is not within the container range
    reference at(geo::TPCID const& id)
      {
        if (hasTPC(id)) return operator[](id);
        throw std::out_of_range("No data for " + std::string(id));
      }

    /// Returns the element for the specified TPC (read-only).
    /// @throw std::out_of_range if the TPC id is not within the container range
    const_reference at (geo::TPCID const& id) const
      {
        if (hasTPC(id)) return operator[](id);
        throw std::out_of_range("No data for " + std::string(id));
      }


    /// Returns the element for the first TPC (unchecked).
    reference first() { return operator[](firstID()); }

    /// Returns the element for the first TPC (unchecked).
    const_reference first() const { return operator[](firstID()); }


    /// Returns the element for the last TPC (unchecked).
    reference last() { return operator[](lastID()); }

    /// Returns the element for the last TPC (unchecked).
    const_reference last() const { return operator[](lastID()); }

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
    void fill(value_type value) { fData.fill(value); }
    
    /// Sets all the elements to a default-constructed `value_type`.
    void clear() { fData.clear(); }
    
    
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
      { return fData.apply(std::forward<Op>(op)); }
    
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
    decltype(auto) apply(Op&& op) const
      { return fData.apply(std::forward<Op>(op)); }
    
  
    /// @}
    // --- END Element access --------------------------------------------------


      private:
    using CryostatNo_t = geo::CryostatID::CryostatID_t;
    using TPCNo_t = geo::TPCID::TPCID_t;

    CryostatNo_t fNCryo; ///< Number of cryostats.
    TPCNo_t fNTPCs; ///< Number of TPCs.
    
    Container_t fData; ///< Data storage.


    /// Returns the internal index of the specified TPC in the storage area.
    size_type index(CryostatNo_t cryo, TPCNo_t tpc) const
      { return (fNTPCs * cryo) + tpc; }

    /// Returns the internal index of the specified TPC in the storage area.
    size_type index(geo::TPCID const& id) const
      { return index(id.Cryostat, id.TPC); }

    /// Returns the ID of the TPC at the specified index (unchecked!)
    geo::TPCID ID(size_type index) const
      { return { index / fNCryo, index % fNCryo }; }

    /// Computes the expected size of this container.
    size_type computeSize() const { return computeSize(fNCryo, fNTPCs); }

    /// Returns the size of a container with the specified dimensions.
    static size_type computeSize(CryostatNo_t nCryo, TPCNo_t nTPCs)
      { return nCryo * nTPCs; }

  }; // class TPCDataContainer<>



  /** **************************************************************************
   * @brief Container with one element per geometry wire plane.
   * @tparam T type of the contained datum
   * @see `geo::GeometryCore::makePlaneData`
   *
   * The container is of fixed size and can't be neither resized nor freed
   * before destruction.
   *
   * This example creates a "map" of tracks starting on each TPC:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto const* geom = lar::providerFrom<geo::GeometryCore>();
   * geo::PlaneDataContainer<std::vector<recob::Hit const*>> hitsPerPlane
   *   (geom->NCryostats(), geom->MaxTPCs(), geom->MaxPlanes());
   *
   * for (recob::Hit const& hit: hits) {
   *   if (hit.WireID()) hitsPerPlane[hit.WireID().planeID()].push_back(&hit);
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
   * * the container assumes the same number of planes in all TPCs, and of TPCs
   *   on all cryostats. While this is not effectively a necessary condition,
   *   keep in mind that this container has no notion whether a given plane or
   *   TPC actually exists in the geometry or not
   * * at least one element is expected to be present
   *
   */
  template <typename T>
  class PlaneDataContainer {

    /// Type of data container helper.
    using Container_t = details::GeoContainerData<T>;
    

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
     * @param nCryo number of cryostats
     * @param nTPCs number of TPCs per cryostat
     * @param nPlanes number of planes per TPC
     *
     * The container is sized to host data for `nCryo` cryostats, each with
     * `nTPCs` TPCs, each one with `fNPlanes` wire planes. Each element in the
     * container is default-constructed.
     */
    PlaneDataContainer
      (unsigned int nCryo, unsigned int nTPCs, unsigned int nPlanes)
      : fNCryo(nCryo)
      , fNTPCs(nTPCs)
      , fNPlanes(nPlanes)
      , fData(computeSize())
      { assert(!empty()); }

    /**
     * @brief Prepares the container with copies of the specified default value.
     * @param nCryo number of cryostats
     * @param nTPCs number of TPCs
     * @param nPlanes number of planes per TPC
     * @param defValue the value to be replicated
     *
     * The container is sized to host data for `nCryo` cryostats, each with
     * `nTPCs` TPCs, and each of them with `fNPlanes` planes. Each element in
     * the container is a copy of defValue.
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
     * auto const* geom = lar::providerFrom<geo::GeometryCore>();
     * geo::TPCDataContainer<unsigned int> countPerPlane
     *   (geom->NCryostats(), geom->MaxTPCs(), geom->MaxPlanes(), 0U);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    PlaneDataContainer(
      unsigned int nCryo, unsigned int nTPCs, unsigned int nPlanes,
      T const& defValue
      )
      : fNCryo(nCryo)
      , fNTPCs(nTPCs)
      , fNPlanes(nPlanes)
      , fData(computeSize(), defValue)
      { assert(!empty()); }


    // --- BEGIN Container status query ----------------------------------------
    /// @name Container status query
    /// @{

    /// Returns the number of elements in the container.
    size_type size() const { return fData.size(); }

    /// Returns the number of elements the container has memory for.
    size_type capacity() const { return fData.capacity(); }

    /// Returns whether the container has no elements (`false` by assumptions).
    bool empty() const { return fData.empty(); }

    /// Returns whether this container hosts data for the specified cryostat.
    bool hasCryostat(geo::CryostatID const& cryoid) const
      { return Container_t::bounded(cryoid.Cryostat, fNCryo); }

    /// Returns whether this container hosts data for the specified TPC.
    bool hasTPC(geo::TPCID const& id) const
      { return hasCryostat(id) && Container_t::bounded(id.TPC, fNTPCs); }

    /// Returns whether this container hosts data for the specified plane.
    bool hasPlane(geo::PlaneID const& id) const
      { return hasTPC(id) && Container_t::bounded(id.Plane, fNPlanes); }

    /// Returns the ID of the first TPC.
    geo::PlaneID firstID() const { return { 0U, 0U, 0U }; }

    /// Returns the ID of the last covered TPC.
    geo::PlaneID lastID() const
      { return { fNCryo - 1U, fNTPCs - 1U, fNPlanes - 1U }; }

    /// @}
    // --- END Container status query ------------------------------------------


    // --- BEGIN Element access ------------------------------------------------
    /// @name Element access
    /// @{

    /// Returns the element for the specified wire plane.
    reference operator[](geo::PlaneID const& id)
      { return fData[index(id)]; }

    /// Returns the element for the specified wire plane (read-only).
    const_reference operator[](geo::PlaneID const& id) const
      { return fData[index(id)]; }

    /// Returns the element for the specified wire plane.
    /// @throw std::out_of_range if the id is not within the container range
    reference at(geo::PlaneID const& id)
      {
        if (hasPlane(id)) return operator[](id);
        throw std::out_of_range("No data for " + std::string(id));
      }

    /// Returns the element for the specified wire plane (read-only).
    /// @throw std::out_of_range if the id is not within the container range
    const_reference at(geo::PlaneID const& id) const
      {
        if (hasPlane(id)) return operator[](id);
        throw std::out_of_range("No data for " + std::string(id));
      }


    /// Returns the element for the first TPC (unchecked).
    reference first() { return operator[](firstID()); }

    /// Returns the element for the first TPC (unchecked).
    const_reference first() const { return operator[](firstID()); }


    /// Returns the element for the last TPC (unchecked).
    reference last() { return operator[](lastID()); }

    /// Returns the element for the last TPC (unchecked).
    const_reference last() const { return operator[](lastID()); }

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
    void fill(value_type value) { fData.fill(value); }
    
    /// Sets all the elements to a default-constructed `value_type`.
    void clear() { fData.clear(); }
    
    
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
      { return fData.apply(std::forward<Op>(op)); }
    
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
    decltype(auto) apply(Op&& op) const
      { return fData.apply(std::forward<Op>(op)); }
    
    
    // --- END Element access --------------------------------------------------


      private:
    using CryostatNo_t = geo::CryostatID::CryostatID_t;
    using TPCNo_t = geo::TPCID::TPCID_t;
    using PlaneNo_t = geo::PlaneID::PlaneID_t;

    CryostatNo_t fNCryo; ///< Number of cryostats.
    TPCNo_t fNTPCs; ///< Number of TPCs per cryostat.
    PlaneNo_t fNPlanes; ///< Number of planes per TPC.

    Container_t fData; ///< Data storage area.

    /// Returns the internal index of the specified TPC in the storage area.
    size_type index(CryostatNo_t cryo, TPCNo_t tpc, PlaneNo_t plane) const
      { return ((cryo) * fNTPCs + tpc) * fNPlanes + plane; }

    /// Returns the internal index of the specified TPC in the storage area.
    size_type index(geo::PlaneID const& id) const
      { return index(id.Cryostat, id.TPC, id.Plane); }

    /// Returns the ID of the TPC at the specified index (unchecked!)
    geo::PlaneID ID(size_type index) const
      {
        return
          { index / fNCryo, index % fNCryo / fNPlanes, index % fNPlanes };
      }

    /// Computes the expected size of this container.
    size_type computeSize() const
      { return computeSize(fNCryo, fNTPCs, fNPlanes); }

    /// Returns the size of a container with the specified dimensions.
    static size_type computeSize
      (CryostatNo_t nCryo, TPCNo_t nTPCs, PlaneNo_t nPlanes)
      { return nCryo * nTPCs * nPlanes; }

  }; // class PlaneDataContainer


  /// @}
  // --- END Geometry data containers ------------------------------------------

} // namespace geo


//------------------------------------------------------------------------------
//--- Template implementation
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

#endif // LARCOREALG_GEOMETRY_GEOMETRYDATACONTAINERS_H
