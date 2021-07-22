/**
 * @file   LineClosestPoint_test.cc
 * @brief  Test of `LineClosestPoint.h` utilities.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   June 29, 2021
 * @see    `larcorealg/Geometry/LineClosestPoint.h`
 */


// Boost libraries
#define BOOST_TEST_MODULE ( LineClosestPoint_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_CLOSE()

// LArSoft libraries
#include "larcorealg/Geometry/LineClosestPoint.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// C++ standard library
#include <utility> // std::pair<>
#include <type_traits> // std::is_same_v<>
#include <cmath> // std::sqrt()


// =============================================================================
void LineClosestPointSimple_test() {
  
  constexpr double tol = 0.001; // percent in CLOSE, absolute in SMALL!
  
  geo::Point_t const p = geo::LineClosestPoint(
    geo::origin() - geo::Zaxis() - 3.0 * geo::Xaxis(), geo::Xaxis(),
    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(), geo::Yaxis()
    );
  
  BOOST_CHECK_SMALL(p.X(),       tol);
  BOOST_CHECK_SMALL(p.Y(),       tol);
  BOOST_CHECK_CLOSE(p.Z(), -1.0, tol);
  
} // LineClosestPointSimple_test()


// -----------------------------------------------------------------------------
void LineClosestPointSimple45_test() {
  
  constexpr double tol = 0.001; // absolute
  
  geo::Point_t const p = geo::LineClosestPoint(
    geo::origin(),                geo::Xaxis(),
    geo::origin() + geo::Yaxis(), (geo::Xaxis() + geo::Zaxis()) / std::sqrt(2.0)
    );
  
  BOOST_CHECK_SMALL(p.X(), tol);
  BOOST_CHECK_SMALL(p.Y(), tol);
  BOOST_CHECK_SMALL(p.Z(), tol);
  
} // LineClosestPointSimple45_test()


// -----------------------------------------------------------------------------
void LineClosestPointAndOffsets_test() {
  
  constexpr double tol = 0.001; // percent in CLOSE, absolute in SMALL!
  
  auto const [ p, ofsA, ofsB ] = geo::LineClosestPointAndOffsets(
    geo::origin()                - 3.0 * geo::Xaxis(), geo::Xaxis(),
    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(), geo::Yaxis()
    );
  
  BOOST_CHECK_SMALL(p.X(), tol);
  BOOST_CHECK_SMALL(p.Y(), tol);
  BOOST_CHECK_SMALL(p.Z(), tol);
  BOOST_CHECK_CLOSE(ofsA,  +3.0, tol);
  BOOST_CHECK_CLOSE(ofsB, -2.0, tol);
  
} // LineClosestPointAndOffsets_test()


// -----------------------------------------------------------------------------
void LineClosestPointWithScaledDirs_test() {
  
  constexpr double tol = 0.001; // percent in CLOSE, absolute in SMALL!
  
  auto const [ p, ofsA, ofsB ] = geo::LineClosestPointAndOffsets(
    geo::origin()                - 3.0 * geo::Xaxis(),  2.0 * geo::Xaxis(),
    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(), -2.0 * geo::Yaxis()
    );
  
  BOOST_CHECK_SMALL(p.X(), tol);
  BOOST_CHECK_SMALL(p.Y(), tol);
  BOOST_CHECK_SMALL(p.Z(), tol);
  BOOST_CHECK_CLOSE(ofsA, +3.0 /  2.0, tol);
  BOOST_CHECK_CLOSE(ofsB, -2.0 / -2.0, tol);
  
} // LineClosestPointWithScaledDirs_test()


// -----------------------------------------------------------------------------
void LineClosestPointWithNonHomogeneousDirs_test() {
  
  constexpr double tol = 0.001; // percent in CLOSE, absolute in SMALL!
  
  auto const [ p, ofsA, ofsB ] = geo::LineClosestPointAndOffsets(
    geo::origin()                - 3.0 * geo::Xaxis(),  1.5 * geo::Xaxis(),
    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(), -2.0 * geo::Yaxis()
    );
  
  BOOST_CHECK_SMALL(p.X(), tol);
  BOOST_CHECK_SMALL(p.Y(), tol);
  BOOST_CHECK_SMALL(p.Z(), tol);
  BOOST_CHECK_CLOSE(ofsA, +3.0 /  1.5, tol);
  BOOST_CHECK_CLOSE(ofsB, -2.0 / -2.0, tol);
  
} // LineClosestPointWithNonHomogeneousDirs_test()


// -----------------------------------------------------------------------------
void LineClosestPointAndOffsetsDocumentation_test() {
  
  constexpr double tol = 0.001; // percent in CLOSE, absolute in SMALL!
  
  /*
   * The promise:
   * 
   * --- 8< --------------------------------------------------------------------
   * The return value is a triplet, which is most easily unpacked immediately:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto [ point, offsetA, offsetB ] = geo::LineClosestPointAndOffsets(
   *   geo::Point_t{ 2, 0, 1 }, geo::Vector_t{ 0.0,   0.5, 0.0 },
   *   geo::Point_t{ 0, 1, 0 }, geo::Vector_t{ 0.866, 0.0, 0.0 }
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will set `point` to `geo::Point{ 2, 1, 1 }`, `offsetA` to `2` and `offsetB`
   * to `2.309...`.
   * To reassign the variables after they have been defined, instead:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::tie(point, offsetA, offsetB) = geo::LineClosestPointAndOffsets(
   *   geo::Point_t{ 0, 1, 0 }, geo::Vector_t{ 0.866, 0.0, 0.0 },
   *   geo::Point_t{ 2, 0, 1 }, geo::Vector_t{ 0.0,   0.5, 0.0 }
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (`point` to `geo::Point{ 2, 1, 0 }`, `offsetA` to `2.039...` and `offsetB`
   * to `2`, because the intersection point is always on the first line).
   * --- 8< --------------------------------------------------------------------
   */
  
  auto [ point, offsetA, offsetB ] = geo::LineClosestPointAndOffsets(
    geo::Point_t{ 2, 0, 1 }, geo::Vector_t{ 0.0,   0.5, 0.0 },
    geo::Point_t{ 0, 1, 0 }, geo::Vector_t{ 0.866, 0.0, 0.0 }
    );
  
  // a way to check we did not mess with the assignment above too much
  static_assert
    (std::is_same_v<decltype(point), geo::Point_t>, "Unexpected point type");
  static_assert
    (std::is_same_v<decltype(offsetA), double>, "Unexpected first offset type");
  static_assert
    (std::is_same_v<decltype(offsetB), double>, "Unexpected second offset type");
  
  BOOST_CHECK_CLOSE(point.X(), 2.0, tol);
  BOOST_CHECK_CLOSE(point.Y(), 1.0, tol);
  BOOST_CHECK_CLOSE(point.Z(), 1.0, tol);
  BOOST_CHECK_CLOSE(offsetA, 2.0, tol);
  BOOST_CHECK_CLOSE(offsetB, 2.0/0.866, tol);
  
  std::tie(point, offsetA, offsetB) = geo::LineClosestPointAndOffsets(
    geo::Point_t{ 0, 1, 0 }, geo::Vector_t{ 0.866, 0.0, 0.0 },
    geo::Point_t{ 2, 0, 1 }, geo::Vector_t{ 0.0,   0.5, 0.0 }
    );
  
  BOOST_CHECK_CLOSE(point.X(), 2.0, tol);
  BOOST_CHECK_CLOSE(point.Y(), 1.0, tol);
  BOOST_CHECK_CLOSE(point.Z(), 0.0, tol);
  BOOST_CHECK_CLOSE(offsetA, 2.0/0.866, tol);
  BOOST_CHECK_CLOSE(offsetB, 2.0, tol);
  
} // LineClosestPointAndOffsetsDocumentation_test()


// =============================================================================
void LineClosestPointWithUnitVectorsSimple_test() {
  
  constexpr double tol = 0.001; // percent in CLOSE, absolute in SMALL!
  
  geo::Point_t const p = geo::LineClosestPointWithUnitVectors(
    geo::origin() - geo::Zaxis() - 3.0 * geo::Xaxis(), geo::Xaxis(),
    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(), geo::Yaxis()
    );
  
  BOOST_CHECK_SMALL(p.X(),       tol);
  BOOST_CHECK_SMALL(p.Y(),       tol);
  BOOST_CHECK_CLOSE(p.Z(), -1.0, tol);
  
} // LineClosestPointWithUnitVectorsSimple_test()


// -----------------------------------------------------------------------------
void LineClosestPointWithUnitVectorsSimple45_test() {
  
  constexpr double tol = 0.001; // absolute
  
  geo::Point_t const p = geo::LineClosestPointWithUnitVectors(
    geo::origin(),                geo::Xaxis(),
    geo::origin() + geo::Yaxis(), (geo::Xaxis() + geo::Zaxis()) / std::sqrt(2.0)
    );
  
  BOOST_CHECK_SMALL(p.X(), tol);
  BOOST_CHECK_SMALL(p.Y(), tol);
  BOOST_CHECK_SMALL(p.Z(), tol);
  
} // LineClosestPointWithUnitVectorsSimple45_test()


// -----------------------------------------------------------------------------
void LineClosestPointWithUnitVectorsAndOffsets_test() {
  
  constexpr double tol = 0.001; // percent in CLOSE, absolute in SMALL!
  
  auto const [ p, ofsA, ofsB ] = geo::LineClosestPointAndOffsetsWithUnitVectors(
    geo::origin()                - 3.0 * geo::Xaxis(), geo::Xaxis(),
    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(), geo::Yaxis()
    );
  
  BOOST_CHECK_SMALL(p.X(), tol);
  BOOST_CHECK_SMALL(p.Y(), tol);
  BOOST_CHECK_SMALL(p.Z(), tol);
  BOOST_CHECK_CLOSE(ofsA, +3.0, tol);
  BOOST_CHECK_CLOSE(ofsB, -2.0, tol);
  
} // LineClosestPointWithUnitVectorsAndOffsets_test()


// -----------------------------------------------------------------------------
void LineClosestPointAndOffsetsWithUnitVectorsDocumentation_test() {
  
  constexpr double tol = 0.001; // percent in CLOSE, absolute in SMALL!
  
  /*
   * The promise:
   * 
   * --- 8< --------------------------------------------------------------------
   * The return value is a triplet, which is most easily unpacked immediately:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * auto [ point, offsetA, offsetB ] = geo::LineClosestPointAndOffsetsWithUnitVectors(
   *   geo::Point_t{ 2, 0, 1 }, geo::Vector_t{ 0, 1, 0 },
   *   geo::Point_t{ 0, 1, 0 }, geo::Vector_t{ 1, 0, 0 }
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * will set `point` to `geo::Point{ 2, 1, 1 }`, `offsetA` to `1` and `offsetB`
   * to `2`.
   * To reassign the variables after they have been defined, instead:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::tie(point, offsetA, offsetB) = geo::LineClosestPointAndOffsetsWithUnitVectors(
   *   geo::Point_t{ 0, 1, 0 }, geo::Vector_t{ 1, 0, 0 },
   *   geo::Point_t{ 2, 0, 1 }, geo::Vector_t{ 0, 1, 0 }
   *   );
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (`point` to `geo::Point{ 2, 1, 0 }`, `offsetA` to `2` and `offsetB` to `1`,
   * because the intersection point is always on the first line).
   * --- 8< --------------------------------------------------------------------
   */
  
  auto [ point, offsetA, offsetB ] = geo::LineClosestPointAndOffsetsWithUnitVectors(
    geo::Point_t{ 2, 0, 1 }, geo::Vector_t{ 0, 1, 0 },
    geo::Point_t{ 0, 1, 0 }, geo::Vector_t{ 1, 0, 0 }
    );
  
  // a way to check we did not mess with the assignment above too much
  static_assert
    (std::is_same_v<decltype(point), geo::Point_t>, "Unexpected point type");
  static_assert
    (std::is_same_v<decltype(offsetA), double>, "Unexpected first offset type");
  static_assert
    (std::is_same_v<decltype(offsetB), double>, "Unexpected second offset type");
  
  BOOST_CHECK_CLOSE(point.X(), 2.0, tol);
  BOOST_CHECK_CLOSE(point.Y(), 1.0, tol);
  BOOST_CHECK_CLOSE(point.Z(), 1.0, tol);
  BOOST_CHECK_CLOSE(offsetA, 1.0, tol);
  BOOST_CHECK_CLOSE(offsetB, 2.0, tol);
  
  std::tie(point, offsetA, offsetB) = geo::LineClosestPointAndOffsetsWithUnitVectors(
    geo::Point_t{ 0, 1, 0 }, geo::Vector_t{ 1, 0, 0 },
    geo::Point_t{ 2, 0, 1 }, geo::Vector_t{ 0, 1, 0 }
    );
  
  BOOST_CHECK_CLOSE(point.X(), 2.0, tol);
  BOOST_CHECK_CLOSE(point.Y(), 1.0, tol);
  BOOST_CHECK_CLOSE(point.Z(), 0.0, tol);
  BOOST_CHECK_CLOSE(offsetA, 2.0, tol);
  BOOST_CHECK_CLOSE(offsetB, 1.0, tol);
  
} // LineClosestPointAndOffsetsWithUnitVectorsDocumentation_test()


// =============================================================================
BOOST_AUTO_TEST_CASE(LineClosestPointTestCase) {
  
  LineClosestPointSimple_test();
  LineClosestPointSimple45_test();
  LineClosestPointAndOffsets_test();
  LineClosestPointWithScaledDirs_test();
  LineClosestPointWithNonHomogeneousDirs_test();
  LineClosestPointAndOffsetsDocumentation_test();
  
} // BOOST_AUTO_TEST_CASE(LineClosestPointTestCase)


// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(LineClosestPointWithUnitVectorsTestCase) {
  
  LineClosestPointWithUnitVectorsSimple_test();
  LineClosestPointWithUnitVectorsSimple45_test();
  LineClosestPointWithUnitVectorsAndOffsets_test();
  LineClosestPointAndOffsetsWithUnitVectorsDocumentation_test();
  
} // BOOST_AUTO_TEST_CASE(LineClosestPointWithUnitVectorsTestCase)


// -----------------------------------------------------------------------------
