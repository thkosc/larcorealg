/**
 * @file   driftvolumes_test.cc
 * @brief  Test for neighbourhood discovery in simple geometries.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 13, 2017
 * 
 * Usage:
 *     
 *     driftvolumes_test configuration.fcl
 *     
 * The configuration file must contain a geometry service provider configuration
 * under `services.Geometry`, That geometry configuration must use the standard
 * channel mapping.
 */

// LArSoft libraries
#include "test/Geometry/geometry_unit_test_base.h"
#include "larcorealg/Geometry/ChannelMapStandardAlg.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcorealg/Geometry/DriftPartitions.h" // BuildDriftVolumes()

// utility libraries
#include "messagefacility/MessageLogger/MessageLogger.h"


#if 0
//------------------------------------------------------------------------------


template <typename T, typename CONT>
void append(std::vector<T>& dest, CONT const& src) {
  using std::cbegin;
  using std::cend;
  dest.insert(dest.end(), cbegin(src), cend(src));
} // append()

template <typename T>
void append(std::vector<T>& dest, std::vector<T>&& src) {
  if (dest.empty()) dest = std::move(src);
  else append(dest, src);
} // append(vector<T>&&)


template <typename Cont>
auto keys(Cont const& m) {
  // TODO use range library instead (this temporary implementation copies keys)
  std::vector<typename Cont::value_type::first_type> result;
  result.reserve(m.size());
  std::transform(
    m.cbegin(), m.cend(),
    std::back_inserter(result), [](auto const& p){ return p.first; }
    );
  return result;

} // keys()

//******************************************************************************
//***  Fundamental classes
//------------------------------------------------------------------------------

template <typename Part>
class PartitionForwardIterator;

template <typename Data>
struct PartitionDataDescriber {
  template <typename Stream>
  PartitionDataDescriber(
    Stream&& out, Data const* data,
    std::string indent = "", std::string firstIndent = ""
    )
    {
      out << firstIndent;
      std::string typeName = lar::debug::demangle<Data>();
      if (data) {
        out << typeName << "[" << ((void*) data) << "]";
      }
      else { 
        out << "no '" << typeName << "' data";
      }
    }
}; // PartitionDataDescriber

/// Describes a `geo::TPCGeo` object for `Partition::describe()`.
template <>
struct PartitionDataDescriber<geo::TPCGeo> {
  template <typename Stream>
  PartitionDataDescriber(
    Stream&& out, geo::TPCGeo const* pTPC,
    std::string indent = "", std::string firstIndent = ""
    )
    {
      out << firstIndent;
      if (!pTPC) {
        out << "no TPC";
        return;
      }
      // verbosity: 2 => ID, size, drift information;
      // 5: also number of wires and active volume
      pTPC->PrintTPCInfo(std::forward<Stream>(out), indent, 2U);
    }
}; // PartitionDataDescriber(TPCGeo)


/**
 * @brief Describes a data object for `Partition::describe()` method.
 * @tparam Stream type of output stream
 * @tparam Data type of data to be described
 * @param out output stream
 * @param data pointer to the data to be described
 * @param indent indentation string applied to all lines except the first
 * @param firstIndent special indentation string for the first line
 * 
 * This function relies on PartitionDataDescriber template class, that should be
 * specialized to customize the printout of known data types.
 */
template <typename Stream, typename Data>
void describePartitionData(
  Stream&& out, Data const* data,
  std::string indent = "", std::string firstIndent = ""
  )
  {
    PartitionDataDescriber<Data>
      (std::forward<Stream>(out), data, indent, firstIndent);
  }



//------------------------------------------------------------------------------
template <typename T>
auto static_less(T a, T b) { return std::less<T>()(a, b); }

template <typename T, typename Key,
  Key KeyExtractor(T const&),
  bool KeyComparer(Key, Key) = static_less<Key>
  >
struct Comparer {
  using Key_t = Key; ///< Type of comparison key.
  using Object_t = T; ///< Type of object to be compared.
  
  bool operator() (Key_t a, Key_t b) const { return key_comp(a, b); }
  bool operator() (Object_t const& a, Key_t b) const
    { return this->operator() (key(a), b); }
  bool operator() (Key_t a, Object_t const& b) const
    { return this->operator() (a, key(b)); }
  bool operator() (Object_t const& a, Object_t const& b) const
    { return this->operator() (key(a), b); }
  
    private:
  static auto key(T const& v) { return KeyExtractor(v); }
  static auto key_comp(Key_t a, Key_t b) { return KeyComparer(a, b); }
  
}; // struct Comparer<>


//------------------------------------------------------------------------------
/// Non-template definitions and data for Partition class hierarchy.
class PartitionBase {
  
    public:
  /// Type of area covered by the partition.
  using Area_t = lar::util::simple_geo::Rectangle<double>;
  
  /// Type of pointer to Area_t data member of type Range_t
  using AreaRangeMember_t
    = PartitionBase::Area_t::Range_t (PartitionBase::Area_t::*);
  
  /// Constructor: sets the covered area and no subpartitions.
  PartitionBase(Area_t const& area): myArea(area) {}
  
  /// Returns whether the specified point is covered by this partition.
  bool contains(double w, double d) const
    { return area().contains(w, d); }
  
  /// Returns the covered area.
  Area_t const& area() const { return myArea; }
  
    protected:
  
  /// Returns a description of the partition area.
  std::string describeArea
    (std::string indent, std::string firstIndent) const;
  
    private:
  Area_t myArea; ///< Covered area.
  
}; // class PartitionBase

std::string PartitionBase::describeArea
  (std::string, std::string firstIndent) const
{
  std::ostringstream sstr;
  sstr << firstIndent << "partition covers " << area();
  return sstr.str();
} // PartitionBase::describeArea()


/// Class extracting the lower bound of the specified range of an area.
template <PartitionBase::AreaRangeMember_t Range>
struct RangeLowerBoundExtractor {
  static constexpr auto range = Range; 
  using Area_t = PartitionBase::Area_t;
  
  double operator() (double lower) const { return lower; }
  double operator() (Area_t::Range_t const& r) const 
    { return (*this)(r.lower); }
  double operator() (Area_t const& area) const 
    { return (*this)(area.*range); }
  
}; // struct RangeLowerBoundExtractor<>


namespace details {
  
  template <typename, typename = void>
  struct is_partition_type: public std::false_type {};
  
  template <typename Part>
  struct is_partition_type
    <
      Part,
      std::enable_if_t
        <std::is_base_of<PartitionBase, std::decay_t<Part>>::value>
    >
    : public std::true_type
  {};
  
  template <typename T>
  constexpr bool is_partition_type_v = is_partition_type<T>();
  
  template <typename, typename = void>
  struct is_partition_ptr: public std::false_type {};
  
  template <typename PartPtr>
  struct is_partition_ptr
    <
      PartPtr,
      std::enable_if_t<is_partition_type_v<decltype(*std::declval<PartPtr>())>>
    >
    : public std::true_type
  {};
  
  template <typename T>
  constexpr bool is_partition_ptr_v = is_partition_ptr<T>();
  
  template <typename, typename = void>
  struct is_partition_ptr_iterator: public std::false_type {};
  
  template <typename Iter>
  struct is_partition_ptr_iterator
    <
      Iter,
      std::enable_if_t
        <is_partition_ptr_v<std::decay_t<typename Iter::value_type>>>
    >
    : public std::true_type
  {};
  
} // namespace details


/// Class extracting the lower bound of the specified range of a partition area.
template <PartitionBase::AreaRangeMember_t Range>
struct PartitionRangeLowerBoundExtractor
  : public RangeLowerBoundExtractor<Range>
{
  using Base_t = RangeLowerBoundExtractor<Range>;
  using Area_t = typename Base_t::Area_t;
  
  using Partition_t = PartitionBase; ///< Base type of the partition.
  
  using Base_t::operator(); // import inherited versions
  
  auto operator()(Partition_t const& part)
    { return Base_t::operator()(part.area()); }
  template <
    typename PartPtr,
    typename = std::enable_if_t<details::is_partition_ptr<PartPtr>::value>
    >
  auto operator()(PartPtr const& part)
    { return operator()(*part); }
  
}; // struct PartitionRangeLowerBoundExtractor<>



/// Base element of a partitioned structure.
template <typename Data>
class Partition: public PartitionBase {
  
    public:
  using Data_t = Data; ///< Type of data stored in the partition.
  using Partition_t = Partition<Data>; ///< This type.
  using Area_t = PartitionBase::Area_t; ///< Type of area.
  
  /// Type of list of subpartitions. It needs to preserve polymorphism.
  using Subpartitions_t = std::vector<std::unique_ptr<Partition_t const>>;
  
  /// Constructor: sets the covered area and no subpartitions.
  Partition(Area_t const& area): PartitionBase(area) {}
  
  /// Destructor (default, virtual).
  virtual ~Partition() = default;
  
  /// Returns the datum directly stored (nullptr if none).
  virtual Data_t* data() const { return nullptr; }
  
  /**
   * @brief Returns the (sub)partition including the specified coordinates.
   * @param w width coordinate of the point
   * @param d depth coordinate of the point
   * @return a pointer to the partition, or `nullptr` if none
   * 
   * This method returns a pointer to the datum associated to the area the
   * specified point (`w`, `d`) belongs to.
   * The area is searched for among all the included partitions.
   */
  virtual Data_t* atPoint(double w, double d) const = 0;
  
  /// Returns a description of the partition.
  virtual std::string describe
    (std::string indent, std::string firstIndent) const
    { return PartitionBase::describeArea(indent, firstIndent); }
  
  /// Returns a description of the partition.
  virtual std::string describe(std::string indent = "") const
    { return describe(indent, indent); }
  
  /// Applies `pred` to this partition first, and then to all subpartitions.
  template <typename Pred>
  void walk(Pred&& pred) const { walk(this, pred); }
  
  /// Returns the number of subparts in the partition (0 if a simple element).
  std::size_t nParts() const { return parts().size(); }
  
    protected:
  
  static Subpartitions_t const NoSubparts; ///< Subpartitions (if any).
  
  /// Returns a list of all subpartitions.
  virtual Subpartitions_t const& parts() const { return NoSubparts; }
  
  /// Applies `pred` to start partition first, and then to all subpartitions.
  template <typename Pred>
  static void walk(Partition_t const* start, Pred&& pred);
  
}; // class Partition<>

template <typename Data>
typename Partition<Data>::Subpartitions_t const Partition<Data>::NoSubparts;

template <typename Data>
template <typename Pred>
void Partition<Data>::walk(Partition_t const* start, Pred&& pred) {
  if (!start) return;
  pred(*start);
  
  // recursive implementation
  for (auto const& subPart: start->parts()) subPart->walk(pred);
  
} // Partition<Data>::walk()


//------------------------------------------------------------------------------
/// Set of drift volumes.
class DriftPartitions {
  
    public:
  /// Type of TPC collection for the partition of a single drift volume.
  using TPCPartition_t = Partition<geo::TPCGeo const>;
  
  /// Type for description of drift range.
  using Range_t = lar::util::simple_geo::Range<double>;
  
  /// Data associated to a single drift volume.
  struct Partition_t {
    /// A partition in width and depth.
    std::unique_ptr<TPCPartition_t> partition;
    /// Interval of drift direction covered by the partition.
    Range_t driftCoverage;
    
    Partition_t(std::unique_ptr<TPCPartition_t>&& part, Range_t const& cover)
      : partition(std::move(part)), driftCoverage(cover) {}
    
    bool coversDrift(double drift) const
      { return driftCoverage.contains(drift); }
    
    static double Position(Partition_t const& part)
      { return part.driftCoverage.lower; }
    using Comparer_t = Comparer<Partition_t, double, Partition_t::Position>;
    
  }; // DriftPartitions::Partition_t
  
  using Position_t = geo::vect::Point_t;
  using Direction_t
    = std::decay_t<decltype(std::declval<geo::TPCGeo>().DriftDir())>;
  using Projection_t
    = ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<double>>;
  using DriftDir_t = Direction_t;

  using Decomposer_t = geo::Decomposer<Direction_t, Position_t, Projection_t>;

  std::vector<Partition_t> partitions; ///< All partitions, sorted by position.
  Decomposer_t decomposer; ///< Decomposition on drift, width and depth axes.
  
  /// Constructor: no partition, but sets the main "drift" direction.
  explicit DriftPartitions(Decomposer_t const& decomp): decomposer(decomp) {}
  
  /// Returns which partition contains the specified position (nullptr if none).
  Partition_t const* partitionAt(Position_t const& pos) const
    { return partitionAt(driftCoord(pos)); }
  
  /// Returns which partition contains the specified drift (nullptr if none).
  Partition_t const* partitionAt(double drift) const;
  
  /// Returns which TPC contains the specified position (nullptr if none).
  geo::TPCGeo const* TPCat(Position_t const& pos) const;
  
  /// Adds the specified partition.
  void addPartition(std::unique_ptr<TPCPartition_t>&& part);
  
  /// Printout of the partition information.
  template <typename Stream>
  void print(Stream&& out) const;
  
  /// Returns drift coordinate (in the drift partition-specific frame) of `pos`.
  double driftCoord(Position_t const& pos) const
    { return decomposer.PointNormalComponent(pos); }
  
  /// Computes the coverage of the specified partition in the drift direction.
  Range_t computeCoverage(TPCPartition_t const& TPCpart) const;
  
  
    private:
  auto partitionAfter(double pos)
    {
      return std::upper_bound
        (partitions.begin(), partitions.end(), pos, Partition_t::Comparer_t());
    }
  auto partitionAfter(double pos) const
    {
      return std::upper_bound
        (partitions.begin(), partitions.end(), pos, Partition_t::Comparer_t());
    }
}; // class DriftPartitions


DriftPartitions::Partition_t const* DriftPartitions::partitionAt
  (double drift) const
{
  auto iPart = partitionAfter(drift); // points to partition after the good one
  if (iPart == partitions.cbegin()) return nullptr;
  if (!(--iPart)->coversDrift(drift)) return nullptr; // maybe not that good?
  return &*iPart;
} // DriftPartitions::partitionAt()


geo::TPCGeo const* DriftPartitions::TPCat(Position_t const& pos) const {
  auto comp = decomposer.DecomposePoint(pos);
  auto part = partitionAt(comp.distance);
  return (part && part->partition)
    ? part->partition->atPoint(comp.projection.X(), comp.projection.Y())
    : nullptr;
} // DriftPartitions::TPCat()


void DriftPartitions::addPartition(std::unique_ptr<TPCPartition_t>&& part) {
  auto range = computeCoverage(*part);
  auto iPart = partitionAfter(range.lower);
  partitions.emplace(iPart, std::move(part), range);
} // DriftPartitions::addPartition()


template <typename Stream>
void DriftPartitions::print(Stream&& out) const {
  
  out << partitions.size() << " drift volume partitions:";
  for (auto const& TPCpart: partitions) {
    out << "\n[" << TPCpart.driftCoverage.lower
      << " -- " << TPCpart.driftCoverage.upper << "]: "
      << TPCpart.partition->describe("  ", "");
  } // for
  out << "\n";
} // DriftPartitions::print()


auto DriftPartitions::computeCoverage
  (TPCPartition_t const& TPCpart) const -> Range_t
{
  /*
   * Computes the range of drift covered by the TPCs in the specified partition.
   * The range is currently defined to include any drift distance covered by
   * any single TPC.
   * The drift distance is computed in "absolute" coordinates, meaning that the
   * origin of the drift is at the origin of the global coordinate system.
   * The drift direction is the one stored in this object.
   */
  
  struct CoverageExtractor {
    Range_t coverage;
    
    CoverageExtractor(DriftPartitions::Decomposer_t const& decomp)
      : decomp(decomp) {}
    
    void operator() (TPCPartition_t const& TPCpart)
      {
        geo::TPCGeo const* TPC = TPCpart.data();
        if (!TPC) return;
        includePoint(TPC->GetCathodeCenter());
        includePoint(TPC->LastPlane().GetCenter());
      }
    
      private:
    DriftPartitions::Decomposer_t const& decomp;
    
    double driftCoord(Position_t const& pos) const
      { return decomp.PointNormalComponent(pos); }
    void includePoint(Position_t const& pos)
      { coverage.extendToInclude(driftCoord(pos)); }
  }; // struct CoverageExtractor
  
  CoverageExtractor extractor(decomposer);
  TPCpart.walk(extractor);
  return extractor.coverage;
  
} // DriftPartitions::computeCoverage()


//******************************************************************************
//***  Partition class hierarchy
//------------------------------------------------------------------------------
/// Partition also containing data directly.
template <typename Data>
class PartitionWithData: public Partition<Data> {
  
    public:
  using Base_t = Partition<Data>; ///< Base class.
  using Partition_t = Partition<Data>; ///< Base type of the partition.
  using Data_t = typename Partition_t::Data_t; ///< Type of contained data.
  using Area_t = typename Partition_t::Area_t; ///< Type of covered area.
  
  /// Constructor: sets the covered area and the contained datum
  PartitionWithData(Area_t const& area, Data_t* myData)
    : Base_t(area), myData(myData)
    {}
  
  /// Returns the datum directly stored (nullptr if none).
  virtual Data_t* data() const override { return myData; }
  
  /// Returns the stored datum only if point is covered, `nullptr` otherwise.
  virtual Data_t* atPoint(double w, double d) const override
    { return Base_t::contains(w, d)? myData: nullptr; }
  
  /// Returns a description of the partition.
  virtual std::string describe
    (std::string indent, std::string firstIndent) const override;
  
    private:
  Data_t* myData; ///< The contained datum.
  
}; // class PartitionWithData


template <typename Data>
std::string PartitionWithData<Data>::describe
  (std::string indent, std::string firstIndent) const
{
  std::string msg = Base_t::describe(indent, firstIndent);
  if (data()) {
    std::ostringstream sstr; 
    sstr << ": ";
    describePartitionData(sstr, data(), indent);
    msg += sstr.str();
  }
  else {
    msg += " (no data)";
  }
  return msg;
} // PartitionWithData<Data>::describe()


//------------------------------------------------------------------------------
/// Unpartitioned element ("leaf") of a partitioned area.
template <typename Data>
class PartitionElement: public PartitionWithData<Data> {
  
    public:
  using Base_t = PartitionWithData<Data>; ///< Base class.
  
  // Import constructors
  using Base_t::Base_t;
  
}; // class PartitionElement


//------------------------------------------------------------------------------
/// Partition divided in subparts.
template <typename Data>
class PartitionContainer: public PartitionWithData<Data> {
  
    public:
  using Base_t = PartitionWithData<Data>; ///< Base class.
  using Partition_t = Partition<Data>; ///< Base type of the partition.
  using Data_t = typename Partition_t::Data_t; ///< Type of contained data.
  using Area_t = typename Partition_t::Area_t; ///< Type of covered area.
  using Subpartitions_t = typename Partition_t::Subpartitions_t;
  
  /// Returns the stored datum only if point is covered, `nullptr` otherwise.
  virtual Data_t* atPoint(double w, double d) const override;
  
  /// Describes each of the subpartitions.
  virtual std::string describe
    (std::string indent, std::string firstIndent) const override;
  
  /// Returns the number of contained subpartitions.
  std::size_t size() const { return parts().size(); }
  
    protected:
  Subpartitions_t myParts; ///< List of subpartitions.
  
  /// Returns a list of the subpartitions owned.
  virtual Subpartitions_t const& parts() const override { return myParts; }
  
  /**
   * @brief Constructor: sets the partition.
   * @param area overall area covered
   * @param subpartitions list of subpartitions (pointers)
   * @param defData datum to be returned for points not covered by subpartitions
   * 
   * The subpartitions will be moved from the argument and will be sorted.
   * Note that this will invalidate pointers to the sub-partitions.
   * 
   * It is required and assumed that the subpartitions do not overlap and that
   * the points covered by them are a subset of `area`.
   * Neither of theses requirements is checked.
   */
  PartitionContainer(
    Area_t const& area,
    Subpartitions_t&& subpartitions,
    Data_t* defData = nullptr
    )
    : Base_t(area, defData), myParts(std::move(subpartitions))
    {}
  
  /// Returns the only partition which could contain the specified width.
  virtual Partition_t const* findPart(double w, double d) const = 0;
  
  /// Introduction to the description of the subpartitions.
  virtual std::string describeIntro() const;
  
}; // class PartitionContainer<>


template <typename Data>
auto PartitionContainer<Data>::atPoint(double w, double d) const -> Data_t* {
  if (!Base_t::contains(w, d)) return nullptr; // not our point at all
  // it's ours; see if it belongs to a subpart
  auto part = findPart(w, d);
  return part? part->atPoint(w, d): Base_t::data();
} // PartitionContainer<Data>::atPoint()


template <typename Data>
std::string PartitionContainer<Data>::describe
  (std::string indent, std::string firstIndent) const
{
  std::string msg = firstIndent + describeIntro();
  if (Base_t::data()) {
    std::ostringstream sstr; 
    sstr << ", and ";
    describePartitionData(sstr, Base_t::data(), indent, "");
    msg += sstr.str();
  }
  
  for (auto const& part: parts()) {
    msg += "\n" + indent + " * ";
    msg += part->describe(indent + "  ", "");
  }
  
  return msg;
} // PartitionContainer<Data>::describe()


template <typename Data>
std::string PartitionContainer<Data>::describeIntro() const {
  return std::to_string(parts().size()) + " subpartitions";
} // PartitionContainer<Data>::describeIntro()


//------------------------------------------------------------------------------
/// Partition of area along the depth dimension.
template <typename Data, typename Sorter>
class SortedPartition: public PartitionContainer<Data> {
  
    public:
  using Base_t = PartitionContainer<Data>; ///< Base class.
  using Partition_t = Partition<Data>; ///< Base type of the partition.
  using Sorter_t = Sorter; ///< Type of sorter being used.
  /// Pointer to constant subpartition.
  using PartitionCPtr_t = std::unique_ptr<Partition_t const>;
  using Data_t = typename Partition_t::Data_t; ///< Type of contained data.
  using Area_t = typename Partition_t::Area_t; ///< Type of covered area.
  
  /// Type of list of subpartitions.
  using Subpartitions_t = typename Base_t::Subpartitions_t;
  
  /**
   * @brief Constructor: sets the partition.
   * @param area overall area covered
   * @param subpartitions list of subpartitions (pointers)
   * @param defData datum to be returned for points not covered by subpartitions
   * @param sorter instance of the sorter to be used
   * 
   * The subpartitions will be moved from the argument and will be sorted.
   * Note that this will invalidate pointers to the sub-partitions.
   * 
   * It is required and assumed that the subpartitions do not overlap and that
   * the points covered by them are a subset of `area`.
   * Neither of theses requirements is checked.
   */
  SortedPartition(
    Area_t const& area,
    Subpartitions_t&& subpartitions,
    Data_t* defData = nullptr,
    Sorter_t sorter = {}
    )
    : Base_t(area, std::move(subpartitions), defData), sorter(sorter)
    { initParts(); }
  
    protected:
  
  Sorter_t sorter;
  
  /// Returns the only partition which could contain the specified key value.
  Partition_t const* findPartWithKey(double key) const;
  
  /// Performs initialization on the specified subpartition list.
  void initParts();
  
}; // class SortedPartition<>


template <typename Data, typename Sorter>
auto SortedPartition<Data, Sorter>::findPartWithKey(double key) const
  -> Partition_t const*
{
  auto pbegin = Base_t::parts().cbegin();
  auto iPart = std::upper_bound(pbegin, Base_t::parts().cend(), key, sorter);
  return (iPart == pbegin)? nullptr: (--iPart)->get();
} // SortedPartition<>::findPartWithKey()


template <typename Data, typename Sorter>
void SortedPartition<Data, Sorter>::initParts() {
  /*
   * Initialization tasks:
   * - ensure that the parts are sorted by increasing depth
   * 
   */
  std::sort(Base_t::myParts.begin(), Base_t::myParts.end(), sorter);
  
} // SortedPartition<>::initParts()


//------------------------------------------------------------------------------
/// Ordering class to sort partition by specified range (lower boundary).
template <PartitionBase::AreaRangeMember_t Range>
struct PartitionSorterByAreaRangeLower {
  using Sorter_t = PartitionSorterByAreaRangeLower<Range>;
  
  using KeyExtractor_t = PartitionRangeLowerBoundExtractor<Range>;
  
  static constexpr auto range = Range;
  
  template <typename T>
  static auto key(T const& obj) { return KeyExtractor_t()(obj); }
  
  /// Type of sorting key. In short: `double`.
  using Key_t = decltype(Sorter_t::key(std::declval<PartitionBase>()));
  
  static Key_t key(Key_t k) { return k; } // shortcut
  static bool sortKey(Key_t a, Key_t b) { return a < b; }
  
  template <typename A, typename B>
  bool operator() (A const& a, B const& b) const
    { return sortKey(key(a), key(b)); }
  
}; // PartitionSorterByAreaRangeLower


/// Partition of area along a area range dimension (width or depth).
template<typename Data, PartitionBase::AreaRangeMember_t Range>
class PartitionSortedByRange
  : public SortedPartition<Data, PartitionSorterByAreaRangeLower<Range>>
{
    public:
  /// Base class.
  using Base_t
    = SortedPartition<Data, PartitionSorterByAreaRangeLower<Range>>;
  using Base_t::Base_t; // import inherited constructors
  
}; // class PartitionSortedByRange


/// Partition of area along the depth dimension.
template <typename Data>
class DepthPartition
  : public PartitionSortedByRange<Data, &Partition<Data>::Area_t::depth>
{
    public:
  /// Base class.
  using Base_t = PartitionSortedByRange<Data, &Partition<Data>::Area_t::depth>;
  
  using Base_t::Base_t; // import inherited constructors
  
  /// Returns the only partition which could contain the specified depth.
  virtual typename Base_t::Partition_t const* findPart
    (double /* w */, double d) const override
    { return Base_t::findPartWithKey(d); }
  
    private:
  
  virtual std::string describeIntro() const override;
  
}; // class DepthPartition


template <typename Data>
std::string DepthPartition<Data>::describeIntro() const {
  std::ostringstream sstr;
  sstr
    << Base_t::size() << " partitions along depth covering " << Base_t::area();
  return sstr.str();
} // DepthPartition<Data>::describeIntro()


/// Partition of area along the width dimension.
template <typename Data>
class WidthPartition
  : public PartitionSortedByRange<Data, &Partition<Data>::Area_t::width>
{
    public:
  /// Base class.
  using Base_t = PartitionSortedByRange<Data, &Partition<Data>::Area_t::width>;
  
  using Base_t::Base_t; // inherited constructors
  
  /// Returns the only partition which could contain the specified depth.
  virtual typename Base_t::Partition_t const* findPart
    (double w, double /* d */) const override
    { return Base_t::findPartWithKey(w); }
  
    private:
  
  virtual std::string describeIntro() const override;
  
}; // class WidthPartition


template <typename Data>
std::string WidthPartition<Data>::describeIntro() const {
  std::ostringstream sstr;
  sstr << Base_t::size() << " partitions along width covering "
    << Base_t::area();
  return sstr.str();
} // WidthPartition<Data>::describeIntro()



//******************************************************************************

using TPCandPos_t = std::pair<double, geo::TPCGeo const*>;

struct TPCgroup_t {
  double pos; ///< Common coordinate of the group.
  std::vector<geo::TPCGeo const*> TPCs;
  
  TPCgroup_t(double pos, std::vector<geo::TPCGeo const*>&& TPCs)
    : pos(pos), TPCs(std::move(TPCs)) {}
  
  static double Position(TPCgroup_t const& tpcg) { return tpcg.pos; }
  using Comparer_t = Comparer<TPCgroup_t, double, TPCgroup_t::Position>;
  
}; // struct TPCgroup_t


std::vector
  <std::pair<DriftPartitions::DriftDir_t, std::vector<geo::TPCGeo const*>>>
groupTPCsByDriftDir(geo::CryostatGeo const& cryo)
{
  
  std::vector
    <std::pair<DriftPartitions::DriftDir_t, std::vector<geo::TPCGeo const*>>>
    result;
  
  lar::util::RealComparisons<double> coordIs(1e-4);
  auto vectorIs = lar::util::makeVector3DComparison(coordIs);
  
  auto const nTPCs = cryo.NTPC();
  for (unsigned int iTPC = 0; iTPC < nTPCs; ++iTPC) {
    geo::TPCGeo const& TPC = cryo.TPC(iTPC);
    
    decltype(auto) driftDir = TPC.DriftDir();
    
    std::size_t iGroup = 0;
    for (; iGroup < result.size(); ++iGroup) {
      if (vectorIs.nonEqual(driftDir, result[iGroup].first)) continue;
      result[iGroup].second.push_back(&TPC);
      break;
    } // for
    
    // if we did not find a group yet, make a new one
    if (iGroup == result.size()) {
      result.emplace_back(
        geo::vect::Rounded01(driftDir, coordIs.threshold),
        std::vector<geo::TPCGeo const*>{ &TPC }
        );
    } // if
    
  } // for
  
  return result;
} // groupTPCsByDriftDir()


std::vector<TPCandPos_t> sortTPCsByDriftCoord(
  std::vector<geo::TPCGeo const*> const& TPCs,
  DriftPartitions::Decomposer_t const& decomp
) {
  //
  // Sorting happens by the drift coordinate of the first wire plane;
  // the absolute value of the drift coordinate value is not relevant nor it is
  // well defined.
  // The result preserves that coordinate for further processing (grouping).
  //
  auto const driftCoord = [&decomp](geo::TPCGeo const& TPC)
    { return decomp.PointNormalComponent(TPC.FirstPlane().GetCenter()); };
  
  std::vector<TPCandPos_t> result;
  result.reserve(TPCs.size());
  std::transform(TPCs.cbegin(), TPCs.cend(), std::back_inserter(result),
    [&driftCoord](geo::TPCGeo const* pTPC)
      { return TPCandPos_t(driftCoord(*pTPC), pTPC); }
    );
  // std::pair sorts by first key first, and second key on par
  // (on par, which may happen often, we don't have means to decide here...)
  std::sort(result.begin(), result.end());
  return result;
} // sortTPCsByDriftCoord()



std::vector<TPCgroup_t> groupByDriftCoordinate
  (std::vector<TPCandPos_t> const& TPCs)
{
  //
  // Produces a list of TPC groups. Within each group, the TPCs have similar
  // drift coordinate.
  // Similar is defined arbitrarily as ten times the plane pitch (of the first
  // planes of the TPC).
  //
  if (TPCs.empty()) return {};
  
  geo::TPCGeo const& firstTPC = *(TPCs.front().second);
  // arbitrary 5 cm if the first TPC has only one plane (pixel readout?);
  // protect against the case where planes have the same position
  // (e.g. dual phase)
  double const groupThickness = 10.0
    * std::min(((firstTPC.Nplanes() > 1)? firstTPC.Plane0Pitch(1): 0.5), 0.1);
  
  auto iFirstTPC = TPCs.cbegin(), tend = TPCs.cend();
  
  std::vector<TPCgroup_t> result;
  while (iFirstTPC != tend) {
    double const posEnd = iFirstTPC->first + groupThickness; // not beyond here
    double sumPos = 0.0;
    std::vector<geo::TPCGeo const*> TPCs;
    auto iEndGroup = iFirstTPC;
    do {
      TPCs.push_back(iEndGroup->second);
      sumPos += iEndGroup->first;
      ++iEndGroup;
    } while ((iEndGroup != tend) && (iEndGroup->first < posEnd));
    
    double const averagePos = sumPos / TPCs.size();
    result.emplace_back(averagePos, std::move(TPCs));
    
    iFirstTPC = iEndGroup;
  } // while (outer)
  
  return result;
} // groupByDriftCoordinate()


unsigned int checkTPCcoords(std::vector<geo::TPCGeo const*> const& TPCs) {
  //
  // verify coordinate system consistency between TPCs:
  //   * need to have the same drift direction
  //   * need to have the same drift coordinate
  //
  
  auto iTPC = TPCs.cbegin(), tend = TPCs.cend();
  if (iTPC == tend) {
    mf::LogProblem("GeometryPartitions")
      << "checkTPCcoords() got an empty partition.";
    return 0;
  }
  
  geo::TPCGeo const& refTPC = **iTPC;
  decltype(auto) refDriftDir = refTPC.DriftDir();
  
  auto driftCoord = [&refDriftDir](geo::TPCGeo const& TPC)
    { return geo::vect::Dot(TPC.FirstPlane().GetCenter(), refDriftDir); };
  
  auto const refDriftPos = driftCoord(refTPC);
  
  lar::util::RealComparisons<double> coordIs(1e-5);
  auto vectorIs = lar::util::makeVector3DComparison(coordIs);
  
  unsigned int nErrors = 0U;
  while (++iTPC != tend) {
    geo::TPCGeo const& TPC = **iTPC;
    
    if (vectorIs.nonEqual(TPC.DriftDir(), refDriftDir)) {
      mf::LogProblem("GeometryPartitions")
        << "Incompatible drift directions between " << TPC.ID()
        << " " << lar::dump::vector3D(TPC.DriftDir()) << " and " << refTPC.ID()
        << " " << lar::dump::vector3D(refTPC.DriftDir());
      ++nErrors;
    }
    auto const driftPos = driftCoord(TPC);
    if (coordIs.nonEqual(driftPos, refDriftPos)) {
      mf::LogProblem("GeometryPartitions")
        << "Incompatible drift coordinate between " << TPC.ID()
        << " (" << driftPos << "( and " << refTPC.ID() << " ("
        << refDriftPos << ")";
      ++nErrors;
    }
  } // while
  return nErrors;
} // checkTPCcoords()


template <typename Range>
DriftPartitions::DriftDir_t detectGlobalDriftDir(Range&& directions) {
  //
  // Returns the first of the specified directions (possibly flipped);
  // throws an exception if any of them is not parallel to that one.
  //
  using std::cbegin;
  using std::cend;
  auto iDir = cbegin(directions);
  auto dend = cend(directions);
  if (!(iDir != dend))
    throw std::runtime_error("detectGlobalDriftDir(): no TPCs provided!");
  
  lar::util::RealComparisons<double> comp(1e-5);
  auto compatibleDir = [comp](auto const& a, auto const& b)
    { return comp.equal(std::abs(geo::vect::Dot(a, b)), +1.0); };
  
  auto const dir = *(iDir++);
  for (; iDir != dend; ++iDir) {
    if (!compatibleDir(dir, *iDir)) {
      throw std::runtime_error(
        "Found drift directions not compatible: "
        + lar::dump::vector3D(dir) + " and " + lar::dump::vector3D(*iDir)
        );
    } // if incompatible
  } // for
  
  // mildly prefer positive directions
  return ((dir.X() <= 0.0) && (dir.Y() <= 0.0) && (dir.Z() <= 0.0))? -dir: dir;
} // detectGlobalDriftDir()



Partition<geo::TPCGeo const>::Area_t TPCarea
  (geo::TPCGeo const& TPC, DriftPartitions::Decomposer_t const& decomposer)
{
  // The area is delimited by TPC bounding box
  
  DriftPartitions::Position_t const lower{ TPC.MinX(), TPC.MinY(), TPC.MinZ() };
  DriftPartitions::Position_t const upper{ TPC.MaxX(), TPC.MaxY(), TPC.MaxZ() };
  
  auto const lowerProj = decomposer.ProjectPointOnPlane(lower);
  auto const upperProj = decomposer.ProjectPointOnPlane(upper);
  
  // we ask to sort the ranges, since the reference base may be flipped
  return {
    { lowerProj.X(), upperProj.X(), true },
    { lowerProj.Y(), upperProj.Y(), true }
    };
  
} // TPCarea()


struct TPCwithArea_t {
  using Area_t = Partition<geo::TPCGeo const*>::Area_t;
  
  Area_t area;
  geo::TPCGeo const* TPC = nullptr;
  
  TPCwithArea_t(Area_t area, geo::TPCGeo const* TPC)
    : area(area), TPC(TPC) {}
  
}; // TPCwithArea_t

/// Class applying a comparison function to keys from the arguments.
template <typename Key, typename ExtractKey, typename Comparer = std::less<Key>>
struct SorterByKey {
  using Key_t = Key;
  
  static bool sortKey(Key_t a, Key_t b) { return Comparer()(a, b); }
  static Key_t key(Key_t k) { return k; }
  template <typename Data>
  static Key_t key(Data const& obj) { return ExtractKey()(obj); }
  
  template <typename A, typename B>
  bool operator() (A const& a, B const& b) const
    { return sortKey(key(a), key(b)); }
  
}; // struct SorterByKey

/// Class extracting the specified range from a Area_t.
template <PartitionBase::AreaRangeMember_t Range>
struct RangeExtractor: public RangeLowerBoundExtractor<Range> {
  using Base_t = RangeLowerBoundExtractor<Range>;
  using Area_t = typename Base_t::Area_t;
  
  using Base_t::operator(); // import inherited versions
  
  double operator() (TPCwithArea_t const& obj) const
    { return Base_t::operator()(obj.area); }
  double operator() (TPCwithArea_t const* ptr) const
    { return operator()(*ptr); }
  
}; // struct RangeExtractor

/// Class sorting any datum by a TPC area range lower boundary.
template<PartitionBase::AreaRangeMember_t Range>
using SortTPCareaByAreaRangeLower = SorterByKey<double, RangeExtractor<Range>>;

using SortTPCwithAreaByWidth
  = SortTPCareaByAreaRangeLower<&TPCwithArea_t::Area_t::width>;
using SortTPCwithAreaByDepth
  = SortTPCareaByAreaRangeLower<&TPCwithArea_t::Area_t::depth>;


std::vector<TPCwithArea_t> addAreaToTPCs(
  std::vector<geo::TPCGeo const*> const& TPCs,
  DriftPartitions::Decomposer_t const& decomposer
) {
  
  std::vector<TPCwithArea_t> result;
  result.reserve(TPCs.size());
  
  for (auto const& TPC: TPCs)
    result.emplace_back(TPCarea(*TPC, decomposer), TPC);
  
  return result;
} // addAreaToTPCs()


auto makeTPCPartitionElement(TPCwithArea_t const& TPCinfo)
{
  return std::make_unique<PartitionElement<geo::TPCGeo const>>
    (TPCinfo.area, TPCinfo.TPC);
}

template <typename BeginIter, typename EndIter>
std::unique_ptr<Partition<geo::TPCGeo const>> makeWidthPartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea);
template <typename BeginIter, typename EndIter>
std::unique_ptr<Partition<geo::TPCGeo const>> makeDepthPartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea);
template <typename BeginIter, typename EndIter>
std::unique_ptr<Partition<geo::TPCGeo const>> makePartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea);


template<
  PartitionBase::AreaRangeMember_t sortingRange,
  typename BeginIter, typename EndIter
  >
std::vector<std::vector<TPCwithArea_t const*>::const_iterator>
groupTPCsByRangeCoord
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea)
{
  /*
   * Groups each TPC with all the following ones overlapping in the selected
   * sorting range.
   * The range of TPCs must be already sorted by lower sorting range coordinate.
   * The result is a list of iterators like the ones in input (in fact, the
   * first is always beginTPCwithArea).
   * The iterators are expected to be valid also after this function has
   * returned (lest the result be unusable).
   */
  
  std::vector<std::vector<TPCwithArea_t const*>::const_iterator> groupStart;
  
  // tolerate 1mm overlaps; this is way forgiving, but apparently there are
  // geometries around (DUNE 35t) with overlaps of that size.
  lar::util::RealComparisons<double> coordIs(0.1);
  
  auto gbegin = beginTPCwithArea;
  while (gbegin != endTPCwithArea) {
    
    groupStart.push_back(gbegin);
    
    //
    // detect the end of this group
    //
    auto range = (*gbegin)->area.*sortingRange;
    auto gend = gbegin;
    while (++gend != endTPCwithArea) {
      //
      // check if the sorting range of this TPC (gend) overlaps the accumulated
      // one; since TPCs are sorted by lower end of the range, gend has that one
      // larger than the accumulated one, and overlap happens only if that lower
      // bound is smaller than the upper bound of the accumulated range;
      // we need to avoid rounding errors: close borders are decided as
      // non-overlapping
      //
      if (coordIs.nonSmaller(((*gend)->area.*sortingRange).lower, range.upper)) break;
      range.extendToInclude((*gend)->area.*sortingRange);
    } // while (inner)
    
    // prepare for the next group
    gbegin = gend;
  } // while (outer)
  
  return groupStart;
  
} // groupTPCsByRangeCoord<>()

/// Returns the TPCs sorted by range coordinate, and the group limits.
template<
  PartitionBase::AreaRangeMember_t sortingRange,
  typename BeginIter, typename EndIter
  >
std::pair<
  std::vector<TPCwithArea_t const*>,
  std::vector<std::vector<TPCwithArea_t const*>::const_iterator>
  >
sortAndGroupTPCsByRangeCoord
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea)
{
  //
  // sort by coordinate; work with pointers for convenience
  //
  std::vector<TPCwithArea_t const*> TPCs(beginTPCwithArea, endTPCwithArea);
  if (TPCs.size() <= 1) return {}; // with only one TPC, refuse to operate
  std::sort
    (TPCs.begin(), TPCs.end(), SortTPCareaByAreaRangeLower<sortingRange>());
  
  //
  // group
  //
  std::vector<std::vector<TPCwithArea_t const*>::const_iterator> TPCgroups
    = groupTPCsByRangeCoord<sortingRange>(TPCs.cbegin(), TPCs.cend());
  assert(!TPCgroups.empty());
  
  return { std::move(TPCs), std::move(TPCgroups) };
  
} // sortAndGroupTPCsByRangeCoord()


template <
  typename BeginIter, typename EndIter, typename TPCendIter,
  typename SubpartMaker
  >
Partition<geo::TPCGeo const>::Subpartitions_t createSubpartitions(
  BeginIter itTPCbegin, EndIter itTPCend, TPCendIter TPCend,
  SubpartMaker subpartMaker
) {
  
  Partition<geo::TPCGeo const>::Subpartitions_t subparts;

  // TPCgroups has an iterator to the starting TPC of each group.
  // Iterators refer to the `TPCs` collection.
  // The end iterator of the group is the starting iterator of the next one;
  // the last group includes all the remaining TPCs and its end iterator is
  // the end iterator of `TPCs`.
  auto igbegin = itTPCbegin;
  while (igbegin != itTPCend) {
    auto const gbegin = *igbegin;
    auto const gend = (++igbegin == itTPCend)? TPCend: *igbegin;
    
    //
    // create a partition from the new group
    //
    if (std::distance(gbegin, gend) == 1) {
      subparts.emplace_back(makeTPCPartitionElement(**gbegin));
    }
    else {
      auto subpart = subpartMaker(gbegin, gend);
      if (!subpart) return {}; // failure!!
      
      subparts.emplace_back(std::move(subpart));
    }
  } // while
  return subparts;
} // createSubpartitions()


template <typename BeginIter, typename EndIter>
TPCwithArea_t::Area_t computeTotalArea(BeginIter TPCbegin, EndIter TPCend) {
  auto iTPC = TPCbegin;
  TPCwithArea_t::Area_t totalArea((*iTPC)->area);
  while (++iTPC != TPCend) totalArea.extendToInclude((*iTPC)->area);
  return totalArea;
} // computeTotalArea()


template <
  typename TPCPartitionResultType,
  PartitionBase::AreaRangeMember_t Range,
  typename BeginIter, typename EndIter,
  typename SubpartMaker
  >
std::unique_ptr<Partition<geo::TPCGeo const>> makeSortedPartition(
  BeginIter beginTPCwithArea, EndIter endTPCwithArea,
  SubpartMaker subpartMaker
) {
  /*
   * TPCs in input are arranged into a partition split on width direction.
   * In case of failure, a null pointer is returned.
   * Do not use this function for a single TPC.
   * The algorithm is as follows:
   * 
   * 1. sort the TPCs by width coordinate
   * 2. for each one, group it with all the following ones overlapping in width
   * 3. for each group with:
   *     .1 more than one element:
   *         .1 let the group be partitioned on depth
   *         .2 if that fails, we fail
   *     .2 just one element: create a partition element
   * 4. the resulting partition is the list of one-TPC "groups" and depth-based
   *    subpartitions of many-TPC groups
   * 5. if the result is a single partition, it must be a depth partition;
   *    then, there is no room for a partition along width, and we declare
   *    failure
   */
  
  //
  // sort by coordinate and group TPCs; work with pointers for convenience
  //
  auto const TPCgroupInfo
    = sortAndGroupTPCsByRangeCoord<Range>(beginTPCwithArea, endTPCwithArea);
  std::vector<TPCwithArea_t const*> const& TPCs = TPCgroupInfo.first;
  std::vector<std::vector<TPCwithArea_t const*>::const_iterator> const&
    TPCgroups = TPCgroupInfo.second;
  
  if (TPCs.empty()) return {}; // failure?
  
  //
  // for each group, create a subpartition
  //
  auto subparts = createSubpartitions
    (TPCgroups.cbegin(), TPCgroups.cend(), TPCs.cend(), subpartMaker);
  
  // if we have grouped everything in a single unit, we have not done any good
  if (subparts.size() == 1) return {};
  
  //
  // compute the total area (it might have been merged in a previous loop...)
  //
  auto totalArea = computeTotalArea(TPCs.cbegin(), TPCs.cend());
  
  //
  // construct and return the final partition
  //
  return std::make_unique<TPCPartitionResultType>
    (totalArea, std::move(subparts));
  
} // makeSortedPartition()


template <typename BeginIter, typename EndIter>
std::unique_ptr<Partition<geo::TPCGeo const>> makeWidthPartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea)
{
  return makeSortedPartition
    <WidthPartition<geo::TPCGeo const>, &TPCwithArea_t::Area_t::width>
    (beginTPCwithArea, endTPCwithArea, &makeDepthPartition<BeginIter, EndIter>);
} // makeWidthPartition()


template <typename BeginIter, typename EndIter>
std::unique_ptr<Partition<geo::TPCGeo const>> makeDepthPartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea)
{
  return makeSortedPartition
    <DepthPartition<geo::TPCGeo const>, &TPCwithArea_t::Area_t::depth>
    (beginTPCwithArea, endTPCwithArea, &makeWidthPartition<BeginIter, EndIter>);
} // makeDepthPartition()


/// A container of partitions organised in a width/depth rectangular grid.
template <typename Data>
class GridPartition: public PartitionContainer<Data> {
    public:
  using Partition_t = Partition<Data>; ///< Base type of the partition.
  using Base_t = PartitionContainer<Data>; ///< Base class.
  using Data_t = typename Partition_t::Data_t; ///< Type of contained data.
  using Area_t = typename Partition_t::Area_t; ///< Type of covered area.
  
  /// Type of list of subpartitions.
  using Subpartitions_t = typename Base_t::Subpartitions_t;
  
  /**
   * @brief Creates a partition with a grid of subpartitions.
   * @param area total area covered by this partition
   * @param subpartitions all subpartitions, row by row
   * @param nDepthPartitions number of partitions on depth direction ("rows")
   * @param nWidthPartitions number of partitions on width direction ("columns")
   * @param defData partition data for areas not covered by subpartitions
   * 
   * The content of the collection of subpartitions is stolen.
   * The subpartitions in that collection are expected to be organized by row:
   * (0;0), (0,1), (0,2)... where the first index spans `nDepthPartitions`
   * values and the second one spans `nWidthPartitions` values.
   * 
   */
  GridPartition(
    Area_t const& area,
    Subpartitions_t&& subpartitions,
    unsigned int nDepthPartitions, unsigned int nWidthPartitions, 
    Data_t* defData = nullptr
    );
  
  /// Constructor: autodetects `nWidthPartitions` from number of subpartitions.
  GridPartition(
    Area_t const& area,
    Subpartitions_t&& subpartitions,
    unsigned int nDepthPartitions, Data_t* defData = nullptr
    )
    : GridPartition(
      area, std::move(subpartitions), nDepthPartitions,
      (nDepthPartitions? subpartitions.size()/nDepthPartitions: 0),
      defData
      )
    {}
  
  /// Prints the information about the partition grid.
  virtual std::string describe
    (std::string indent, std::string firstIndent) const override;
  
  
    private:
  std::size_t nWidthParts; ///< Number of partitions on width direction.
  std::size_t nDepthParts; ///< Number of partitions on depth direction.
  std::vector<double> widthSeps; ///< Separators for width dimension.
  std::vector<double> depthSeps; ///< Separators for depth dimension.
  
  auto part(std::size_t iDepth, std::size_t iWidth) ->decltype(auto)
    { return Base_t::parts()[iDepth * nWidthParts + iWidth]; }
  auto part(std::size_t iDepth, std::size_t iWidth) const ->decltype(auto)
    { return Base_t::parts()[iDepth * nWidthParts + iWidth]; }
  
  /// Returns the only partition which could contain the specified depth.
  virtual Partition_t const* findPart(double w, double d) const override;
  
  std::vector<double> computeWidthSeps() const;
  std::vector<double> computeDepthSeps() const;
  
  template<
    PartitionBase::AreaRangeMember_t Range, typename BeginIter, typename EndIter
    >
  static std::vector<double> detectSeparators(
    BeginIter b, EndIter e,
    std::size_t const nGroups,
    std::size_t const startDelta,
    std::size_t const stride
    );
  
}; // class GridPartition<>


template <typename Data>
GridPartition<Data>::GridPartition(
  Area_t const& area,
  Subpartitions_t&& subpartitions,
    unsigned int nDepthPartitions, unsigned int nWidthPartitions, 
  Data_t* defData /* = nullptr */
  )
  : Base_t(area, std::move(subpartitions), defData)
  , nWidthParts(nWidthPartitions)
  , nDepthParts(nDepthPartitions)
  , widthSeps(computeWidthSeps())
  , depthSeps(computeDepthSeps())
{
  assert(nDepthParts * nWidthParts == Base_t::size());
} // GridPartition<Data>::GridPartition()


template <typename Data>
std::string GridPartition<Data>::describe
  (std::string indent, std::string firstIndent) const
{
  std::ostringstream sstr; 
  sstr << firstIndent << Base_t::describeIntro()
    << " in a (WxD) = " << nWidthParts << " x " << nDepthParts << " grid";
  if (Base_t::data()) {
    sstr << ", and ";
    describePartitionData(sstr, Base_t::data(), indent, "");
  }
  for (std::size_t iDepth = 0; iDepth < nDepthParts; ++iDepth) {
    for (std::size_t iWidth = 0; iWidth < nWidthParts; ++iWidth) {
      sstr << "\n" << indent << " [" << iDepth << "][" << iWidth << "] "
        << part(iDepth, iWidth)->describe(indent + "  ", "");
    } // for width
  } // for depth
  
  return sstr.str();
} // GridPartition<Data>::describe()


template <typename Data>
auto GridPartition<Data>::findPart(double w, double d) const
  -> Partition_t const* 
{
  auto const iWidth = std::upper_bound(widthSeps.cbegin(), widthSeps.cend(), w);
  if (iWidth == widthSeps.cbegin()) return nullptr;
  auto const iDepth = std::upper_bound(depthSeps.cbegin(), depthSeps.cend(), d);
  if (iDepth == depthSeps.cbegin()) return nullptr;
  return part(
    std::distance(depthSeps.cbegin(), iDepth) - 1U,
    std::distance(widthSeps.cbegin(), iWidth) - 1U
    ).get();
} // GridPartition<Data>::findPart()


template <typename Data>
std::vector<double> GridPartition<Data>::computeWidthSeps() const {
  return detectSeparators<&Area_t::width>(
    Base_t::parts().cbegin(), Base_t::parts().cend(),
    nWidthParts, 1U, nWidthParts
    );
} // GridPartition<Data>::computeWidthSeps()

template <typename Data>
std::vector<double> GridPartition<Data>::computeDepthSeps() const {
  return detectSeparators<&Area_t::depth>(
    Base_t::parts().cbegin(), Base_t::parts().cend(),
    nDepthParts, nWidthParts, 1U
    );
} // GridPartition<Data>::computeDepthSeps()

template <typename Data>
template
  <PartitionBase::AreaRangeMember_t Range, typename BeginIter, typename EndIter>
std::vector<double> GridPartition<Data>::detectSeparators(
  BeginIter b, EndIter e,
  std::size_t const nGroups,
  std::size_t const startDelta,
  std::size_t const stride
) {
  /*
   * The iterators are better be random access.
   * The range [b,e[ is considered to be a 2D table stored row after row.
   * This function can operate on rows or columns, given the proper arguments.
   * 
   * The full range is split in nGroups "groups" (e.g. rows or columns).
   * Each group g starts at the element g x startDelta.
   * All members of that group are stride elements far from each other.
   * 
   * Be the table of size (d x w).
   * To process data row by row:
   * - the group is a row
   * - the number of groups is the number of rows: nGroups = d
   * - each group will have (d x w)/nGroups = w elements
   * - the start offset is the number of group times the size of it:
   *   startDelta = w
   * - the elements are contiguous: stride = 1
   * 
   * To process data column by column:
   * - the group is a column
   * - the number of groups is the number of columns: nGroups = w
   * - each group will have (d x w)/nGroups = d elements
   * - the start offset matches the number of group: startDelta = 1
   * - the elements are separated by a full row: stride = w
   * 
   */
  // the separator is on the lower bound of selected range of the partition area
  
  static_assert(details::is_partition_ptr_iterator<BeginIter>(),
    "Begin iterator does not point to a pointer to partition type");
  
  using Range_t = PartitionBase::Area_t::Range_t;
  PartitionRangeLowerBoundExtractor<Range> lowerBound;
  
  std::size_t const nParts = std::distance(b, e);
  std::size_t const nPartsInGroup = nParts / nGroups;
  
  auto const part
    = [b](std::size_t index){ return std::next(b, index)->get(); };
  
  std::vector<double> seps(nGroups);
  for (size_t g = 0; g < nGroups; ++g) {
  
    double& sep = seps[g];
    
    // indices of an element in the previous and next group, respectively
    std::size_t index = g * startDelta;
    sep = lowerBound(part(index));
    
    std::size_t const iend = index + nPartsInGroup * stride;
    while ((index += stride) < iend) {
      double const l = lowerBound(part(index));
      if (sep > l) sep = l;
    } // while (element)
    
  } // while (groups)
  return seps;
} // GridPartition<Data>::detectSeparators()


template <typename BeginIter, typename EndIter>
std::unique_ptr<Partition<geo::TPCGeo const>> makeGridPartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea)
{
  /*
   * Requires at least 4 input TPCs (otherwise, do not use GridPartition).
   * 
   * 1. attempt a partition on width
   *    1. if failed, return failure
   * 2. attempt to partition the first subpartition on depth
   *    1. if failed, return failure
   * 3. extend each depth partition to the other width partitions
   *    1. if the extension of one partition line fails, discard it
   *    2. if no partition line survives, return failure
   * 4. run makePartition() on each of the cells
   * 5. create and return a GridPartition object from the subpartitions so
   *    created
   * 
   * This algorithm could use some factorization...
   */
  using Area_t = TPCwithArea_t::Area_t;
  
  //
  // sort by width coordinate; work with pointers for convenience
  //
  auto const TPCgroupInfo = sortAndGroupTPCsByRangeCoord<&Area_t::width>
    (beginTPCwithArea, endTPCwithArea);
  std::vector<TPCwithArea_t const*> const& TPCs = TPCgroupInfo.first;
  std::vector<std::vector<TPCwithArea_t const*>::const_iterator> const&
    TPCgroups = TPCgroupInfo.second;
  
  if (TPCs.empty()) return {}; // failure?
  // with only one TPC, then makeTPCPartitionElement() should be used instead!
  if (TPCs.size() < 4) return {};
  
  unsigned int const nWidthParts = TPCgroups.size();
  if (nWidthParts <= 1) return {}; // only one group ain't no good
  
  //
  // sort TPCs in the first width partition by depth coordinate
  //
  auto const FirstColGroupInfo
    = sortAndGroupTPCsByRangeCoord<&Area_t::depth>(TPCgroups[0], TPCgroups[1]);
  std::vector<TPCwithArea_t const*> const& FirstColTPCs
    = FirstColGroupInfo.first;
  std::vector<std::vector<TPCwithArea_t const*>::const_iterator> const&
    FirstColGroups = FirstColGroupInfo.second;
  
  if (FirstColTPCs.empty()) return {}; // failure?
  if (FirstColGroups.size() <= 1 ) return {}; // only one row ain't good either
  
  //
  // collect all candidate separation ranges
  //
  // First depth partition has no lower limit, last one has no upper limit
  // (they include all TPCs with depth lower than the upper limit in the first
  // case, all TPCs with depth higher of the lower limit in the last case).
  // Checks need to be done in the gaps between depth partitions.
  // So we start by skipping the first border.
  
  lar::util::RealComparisons<double> coordIs(0.1);
  std::vector<Area_t::Range_t> depthGaps; // candidate gaps
  auto icnext = FirstColGroups.cbegin(), icprev = icnext,
    icend = FirstColGroups.cend();
  while (++icnext != icend) {
    //
    // establish the range of the group in the candidate depth partition
    // from the group in the first width partition
    //
    auto const cprev = *icprev;
    auto const cnext = *icnext;
    
    depthGaps.emplace_back
      ((*cprev)->area.depth.upper, (*cnext)->area.depth.lower);
    
    icprev = icnext;
  } // while
  assert(!depthGaps.empty());
  
  //
  // see that for every other width partition separations hold
  //
  auto igbegin = TPCgroups.cbegin();
  while (++igbegin != TPCgroups.cend()) {
    //
    // prepare the TPC groups within this width partition
    //
    auto igend = std::next(igbegin);
    auto gbegin = *igbegin;
    auto gend = (igend == TPCgroups.cend())? TPCs.cend(): *igend;
    
    auto const ColGroupInfo
      = sortAndGroupTPCsByRangeCoord<&Area_t::depth>(gbegin, gend);
    std::vector<TPCwithArea_t const*> const& ColTPCs = ColGroupInfo.first;
    std::vector<std::vector<TPCwithArea_t const*>::const_iterator> const&
      ColGroups = ColGroupInfo.second;
    
    // failure to partition a single column means total failure
    if (ColTPCs.empty()) return {};
    if (ColGroups.size() <= 1) return {}; // only one row ain't good either
    
    //
    // compute the coverage of each of the depth groups
    //
    std::vector<TPCwithArea_t::Area_t::Range_t> groupDepths(ColGroups.size());
    auto iGDepth = groupDepths.begin();
    for (auto icgstart = ColGroups.cbegin(); icgstart != ColGroups.cend();
      ++icgstart, ++iGDepth)
    {
      auto const icgend = std::next(icgstart);
      auto ictpc = *icgstart;
      auto const ictend = (icgend == ColGroups.cend())? ColTPCs.cend(): *icgend;
      while (ictpc != ictend) 
        iGDepth->extendToInclude((*(ictpc++))->area.depth);
    } // for
    
    //
    // check each of the remaining candidate gaps
    //
    auto iGap = depthGaps.begin();
    while (iGap != depthGaps.end()) {
      Area_t::Range_t& gap = *iGap;
      
      //
      // check that the gap holds
      //
      bool bGoodGap = false;
      // first TPC starting after the gap (even immediately after):
      auto iCGroup = std::lower_bound(
        groupDepths.cbegin(), groupDepths.cend(), gap.upper,
        SortTPCwithAreaByDepth()
        );
      
      // any TPCs before/after this gap?
      if ((iCGroup != groupDepths.begin()) && (iCGroup != groupDepths.end())) {
        Area_t::Range_t const& before = *(std::prev(iCGroup));
        Area_t::Range_t const& after = *iCGroup;
        Area_t::Range_t const TPCgap{ before.upper, after.lower };
        
        // correct the gap
        if (coordIs.strictlySmaller(iGap->lower, TPCgap.lower))
          iGap->lower = TPCgap.lower;
        if (coordIs.strictlyGreater(iGap->upper, TPCgap.upper))
          iGap->upper = TPCgap.upper;
        
        // if nothing is left, gap is gone
        bGoodGap = coordIs.nonSmaller(iGap->upper, iGap->lower);
      } // if TPCs around the gap
      
      //
      // if the gap has been flagged as bad, remove it
      //
      if (bGoodGap) ++iGap;
      else iGap = depthGaps.erase(iGap);
      
    } // while (separation)
    
    if (depthGaps.empty()) return {}; // no surviving gaps means failure
    
  } // while (width partition)
  
  // 
  // turn the gaps into separators
  // 
  std::vector<double> depthSep;
  std::transform(
    depthGaps.cbegin(), depthGaps.cend(), std::back_inserter(depthSep),
    [](auto const& r){ return (r.lower + r.upper) / 2.0; }
    );
  unsigned int const nDepthParts = depthSep.size() + 1;
  
  //
  // fill the groups with TPCs, and create subpartitions from each of them
  //
  Partition<geo::TPCGeo const>::Subpartitions_t subparts
    (nWidthParts * nDepthParts);
  Area_t totalArea;
  
  unsigned int iWidth = 0;
  for (auto igbegin = TPCgroups.cbegin(); igbegin != TPCgroups.cend();
    ++igbegin, ++iWidth
  ) {
    
    // sort TPCs in this group (yes, again; this time we don't group just yet)
    auto igend = std::next(igbegin);
    auto gbegin = *igbegin;
    auto gend = (igend == TPCgroups.cend())? TPCs.cend(): *igend;
    std::vector<TPCwithArea_t const*> ColTPCs(gbegin, gend);
    std::sort(ColTPCs.begin(), ColTPCs.end(),
      SortTPCareaByAreaRangeLower<&Area_t::depth>());
    
    unsigned int iDepth = 0;
    auto cgstart = ColTPCs.cbegin(), TPCend = ColTPCs.cend();
    for (double sep: depthSep) {
      
      //
      // collect all TPCs for this partition
      //
      // the first TPC that starts *after* the separator:
      auto cgend
        = std::upper_bound(cgstart, TPCend, sep, SortTPCwithAreaByDepth());
      // if we cut out TPCs that were included because of some tolerance,
      // recover them now
      while (cgend != cgstart) {
        auto cglast = std::prev(cgend);
        if (coordIs.strictlySmaller((*cglast)->area.depth.lower, sep)) break;
        cgend = cglast;
      } // while
      assert(cgstart != cgend); // separator selection should guarantee this
      
      //
      // create and register the partition
      //
      auto part = makePartition(cgstart, cgend);
      if (!part) return {}; // late failure!
      totalArea.extendToInclude(part->area());
      subparts[iDepth * nWidthParts + iWidth] = std::move(part);
      
      ++iDepth;
      cgstart = cgend;
    } // for all depth separators
    
    //
    // collect all the TPCs after the last separator
    //
    auto part = makePartition(cgstart, TPCend);
    if (!part) return {}; // super-late failure!
    totalArea.extendToInclude(part->area());
    subparts[iDepth * nWidthParts + iWidth] = std::move(part);
    
  } // for all width partitions
  
  return std::make_unique<GridPartition<geo::TPCGeo const>>
    (totalArea, std::move(subparts), nWidthParts, nDepthParts);
  
} // makeGridPartition()


template <typename BeginIter, typename EndIter>
std::unique_ptr<Partition<geo::TPCGeo const>> makePartition
  (BeginIter beginTPCwithArea, EndIter endTPCwithArea)
{
  /*
   * Organizes a list of TPCs in a hierarchical partition.
   * Three main elements are used:
   * - single element partition objects: that's the single TPC end point
   * - TPC groups organised along width
   * - TPC groups organised along depth
   * 
   * The procedure is recursively analysing a set of TPCs:
   * - if the set is actually one TPC only, use a PartitionElement
   * - attempt partitioning o a grid; if fails:
   * - attempt partitioning along width:
   *    * determine overlapping groups: a set of TPCs which share part of the
   *      width range
   *    * recurse on each overlapping group with more than one TPC,
   *      attempting a depth partition
   *    * if that fails, bail out since we don't have code to deal with a layout
   *      with areas overlapping on both directions at the same time
   *    * add the single elements and the overlapping groups to the width
   *      partition
   * - attempt partitioning along height:
   *    * same algorithm as for width
   * - pick the partitioning with less elements
   */
  using value_type = std::remove_reference_t<decltype(*beginTPCwithArea)>;
  static_assert(
    std::is_pointer<value_type>()
    && std::is_same
      <std::decay_t<std::remove_pointer_t<value_type>>, TPCwithArea_t>(),
    "Iterators must point to TPCwithArea_t pointers."
    );
  
  
  auto const size = std::distance(beginTPCwithArea, endTPCwithArea);
  if (size == 1) {
    return makeTPCPartitionElement(**beginTPCwithArea);
  }
  
  auto gPart = makeGridPartition(beginTPCwithArea, endTPCwithArea);
  if (gPart) return gPart;
  
  auto wPart = makeWidthPartition(beginTPCwithArea, endTPCwithArea);
  auto dPart = makeDepthPartition(beginTPCwithArea, endTPCwithArea);
  
  if (wPart) {
    
    if (dPart) { // wPart && dPart
      if (wPart->nParts() < dPart->nParts()) return wPart;
      else                                   return dPart; // slight preference
    }
    else { // wPart && !dPart
      return wPart; // easy choice
    }
    
  }
  else {
    
    if (dPart) { // !wPart && dPart
      return dPart; // easy choice
    }
    else { // !wPart && !dPart
      return {}; // failure!!
    }
    
  }
  
} // makePartition(Iter)


template <typename BeginIter, typename EndIter>
auto makeCPointerVector(BeginIter b, EndIter e) {
  using value_type = typename BeginIter::value_type;
  std::vector<value_type const*> result;
  result.reserve(std::distance(b, e));
  std::transform
    (b, e, std::back_inserter(result), std::addressof<value_type const>);
  return result;
} // makeCPointerVector()

template <typename T>
auto makeCPointerVector(std::vector<T> const& v)
  { return makeCPointerVector(v.cbegin(), v.cend()); }


std::unique_ptr<DriftPartitions::TPCPartition_t> makePartition
  (std::vector<TPCwithArea_t> const& TPCs)
{
  // TODO use range library instead:
//  auto TPCptrs = TPCs | ranges::view::transform(std::addressof);
  auto TPCptrs = makeCPointerVector(TPCs);
  using std::cbegin;
  using std::cend;
  return makePartition(cbegin(TPCptrs), cend(TPCptrs)); 
} // makePartition(coll)


DriftPartitions BuildDriftVolumes(geo::CryostatGeo const& cryo) {
  
  //
  // group TPCs by drift direction
  //
  auto TPCsByDriftDir = groupTPCsByDriftDir(cryo);
  
  //
  // determine the cryostat-wide drift direction (arbitrary but consistent)
  // and the decomposition base (using the same for all drift partitions);
  //
  
  // In practice we use the coordinate system from the first TPC;
  // we still check that all drift directions are compatible,
  // but the result of detection is ignored.
  /* auto globalDriftDir = */ detectGlobalDriftDir(keys(TPCsByDriftDir));
  DriftPartitions::Position_t cryoCenter
    { cryo.CenterX(), cryo.CenterY(), cryo.CenterZ() };
  geo::TPCGeo const& firstTPC = cryo.TPC(0);
  DriftPartitions::Decomposer_t decomposer
    ({ cryoCenter, firstTPC.RefWidthDir(), firstTPC.RefDepthDir() });
  
  //
  // further group TPCs by plane position in drift direction
  //
  std::vector<TPCgroup_t> TPCgroups;
  for (auto const& TPCsOnDriftDir: TPCsByDriftDir) {
    auto TPCs = sortTPCsByDriftCoord(TPCsOnDriftDir.second, decomposer);
    append(TPCgroups, groupByDriftCoordinate(TPCs));
  } // for
  
  //
  // verify coordinate system consistency between TPCs
  //
  for (auto const& TPCgroup: TPCgroups) {
    unsigned int errors = checkTPCcoords(TPCgroup.TPCs);
    if (errors > 0) {
      throw std::runtime_error(
        "TPCs in partition have different drift directions ("
        + std::to_string(errors) + " errors found in "
        + std::to_string(TPCgroup.TPCs.size()) + " TPCs)."
        );
    } // if
  } // for
  
  //
  // partition each group 
  //
  DriftPartitions partitions(decomposer);
  for (auto const& TPCgroup: TPCgroups) {
    auto TPCs = addAreaToTPCs(TPCgroup.TPCs, decomposer);
    auto part = makePartition(TPCs);
    if (!part) {
      cet::exception e("BuildDriftVolumes");
      e << "Failed to construct partition out of " << TPCs.size() << " TPCs:";
      for (auto const& TPCinfo: TPCs) {
        e << "\n at " << TPCinfo.area << " TPC ";
        TPCinfo.TPC->PrintTPCInfo(e, "   ", 5U);
      } // for
      throw e;
    } // if error
    partitions.addPartition(std::move(part));
  } // for
  
  return partitions;
} // BuildDriftVolumes()

#endif // 0

//------------------------------------------------------------------------------
//---  The test environment
//---

// we define here all the configuration that is needed;
// we use an existing class provided for this purpose, since our test
// environment allows us to tailor it at run time.
using StandardGeometryConfiguration
  = testing::BasicGeometryEnvironmentConfiguration<geo::ChannelMapStandardAlg>;

/*
 * GeometryTesterFixture, configured with the object above, is used in a
 * non-Boost-unit-test context.
 * It provides:
 * - `geo::GeometryCore const* Geometry()`
 * - `geo::GeometryCore const* GlobalGeometry()` (static member)
 */
using StandardGeometryTestEnvironment
  = testing::GeometryTesterEnvironment<StandardGeometryConfiguration>;


//------------------------------------------------------------------------------
//---  The tests
//---

/** ****************************************************************************
 * @brief Runs the test
 * @param argc number of arguments in argv
 * @param argv arguments to the function
 * @return number of detected errors (0 on success)
 * @throw cet::exception most of error situations throw
 * 
 * The arguments in argv are:
 * 0. name of the executable ("Geometry_test")
 * 1. path to the FHiCL configuration file
 * 2. FHiCL path to the configuration of the geometry
 *    (default: services.Geometry)
 * 
 */
//------------------------------------------------------------------------------
int main(int argc, char const** argv) {
  
  StandardGeometryConfiguration config("geometry_test");
  config.SetMainTesterParameterSetName("geotest");
  
  //
  // parameter parsing
  //
  int iParam = 0;
  
  // first argument: configuration file (mandatory)
  if (++iParam < argc) config.SetConfigurationPath(argv[iParam]);
  
  // second argument: path of the parameter set for geometry configuration
  // (optional; default: "services.Geometry" from the inherited object)
  if (++iParam < argc) config.SetGeometryParameterSetPath(argv[iParam]);
  
  //
  // testing environment setup
  //
  StandardGeometryTestEnvironment TestEnvironment(config);
  auto const& geom = *(TestEnvironment.Provider<geo::GeometryCore>());
  
  //
  // run the test algorithm
  //
  
  unsigned int nErrors = 0;
  for (auto const& cryo: geom.IterateCryostats()) {
    auto partition = geo::buildDriftVolumes(cryo);
    mf::LogVerbatim("driftvolumes_test")
      << "Partition for cryostat " << cryo.ID() << ":";
    partition.print(mf::LogVerbatim("driftvolumes_test"));
    
    //
    // test that the partition topology is correct
    //
    for (geo::TPCGeo const& TPC: geom.IterateTPCs(cryo.ID())) {
      
      decltype(auto) center = TPC.GetCenter();
      
      auto where = partition.TPCat(center);
      if (!where) {
        mf::LogProblem("driftvolumes_test")
          << "Center of TPC " << TPC.ID() << " " << lar::dump::vector3D(center)
          << " not assigned to any TPC!";
        ++nErrors;
      }
      else if (where->ID() != TPC.ID()) {
        mf::LogProblem log("driftvolumes_test");
        log
          << "Center of TPC " << TPC.ID() << " " << lar::dump::vector3D(center)
          << " assigned to TPC " << where->ID() << ":\n";
        where->PrintTPCInfo(log, "  ", /* verbosity */ 5);
        ++nErrors;
      }
      
      //
      // test ownership of a uniform distribution of points inside the TPC
      //
      constexpr int nXsteps = 5, nYsteps = 5, nZsteps = 5;
      const double xstep = TPC.HalfSizeX() / (nXsteps + 1);
      const double ystep = TPC.HalfSizeY() / (nYsteps + 1);
      const double zstep = TPC.HalfSizeZ() / (nZsteps + 1);
      for (int xs = -nXsteps; xs <= nXsteps; ++xs) {
        double const x = center.X() + xs * xstep;
        for (int ys = -nYsteps; ys <= nYsteps; ++ys) {
          double const y = center.Y() + ys * ystep;
          for (int zs = -nZsteps; zs <= nZsteps; ++zs) {
            double const z = center.Z() + zs * zstep;
            
            auto where = partition.TPCat({ x, y, z});
            if (!where) {
              mf::LogProblem("driftvolumes_test")
                << "Point " << lar::dump::vector3D(center) << " within TPC "
                << TPC.ID() << " (" << xs << "," << ys << "," << zs
                << ") not assigned to any TPC!";
              ++nErrors;
            }
            else if (where->ID() != TPC.ID()) {
              mf::LogProblem log("driftvolumes_test");
              log
                << "Point " << lar::dump::vector3D(center) << " within TPC "
                << TPC.ID() << " (" << xs << "," << ys << "," << zs
                << ") assigned to TPC " << where->ID() << ":\n";
              where->PrintTPCInfo(log, "  ", /* verbosity */ 5);
              ++nErrors;
            }
            
          } // for zs
        } // for ys
      } // for xs
      
    } // for TPCs
    
  } // for
  
  // 4. And finally we cross fingers.
  if (nErrors > 0) {
    mf::LogError("geometry_test") << nErrors << " errors detected!";
  }
  
  return nErrors;
} // main()

