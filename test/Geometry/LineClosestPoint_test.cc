/**
 * @file   LineClosestPoint_test.cc
 * @brief  Test of `LineClosestPoint.h` utilities.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   June 29, 2021
 * @see    `larcorealg/Geometry/LineClosestPoint.h`
 */

// Boost libraries
#define BOOST_TEST_MODULE LineClosestPoint_test
#include <boost/test/unit_test.hpp> // BOOST_AUTO_TEST_CASE(), BOOST_TEST()

// LArSoft libraries
#include "larcorealg/Geometry/LineClosestPoint.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_vectors.h"

// C++ standard library
#include <cmath>       // std::sqrt()
#include <type_traits> // std::is_same_v<>
#include <utility>     // std::pair<>

// =============================================================================
void LineClosestPointSimple_test()
{

  auto const tol = boost::test_tools::tolerance(0.001);

  geo::Point_t const p = geo::LineClosestPoint(geo::origin() - geo::Zaxis() - 3.0 * geo::Xaxis(),
                                               geo::Xaxis(),
                                               geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(),
                                               geo::Yaxis());

  BOOST_TEST(p.X() == 0.0, tol);
  BOOST_TEST(p.Y() == 0.0, tol);
  BOOST_TEST(p.Z() == -1.0, tol);

} // LineClosestPointSimple_test()

// -----------------------------------------------------------------------------
void LineClosestPointSimple45_test()
{

  auto const tol = boost::test_tools::tolerance(0.001);

  geo::Point_t const p = geo::LineClosestPoint(geo::origin(),
                                               geo::Xaxis(),
                                               geo::origin() + geo::Yaxis(),
                                               (geo::Xaxis() + geo::Zaxis()) / std::sqrt(2.0));

  BOOST_TEST(p.X() == 0.0, tol);
  BOOST_TEST(p.Y() == 0.0, tol);
  BOOST_TEST(p.Z() == 0.0, tol);

} // LineClosestPointSimple45_test()

// -----------------------------------------------------------------------------
void LineClosestPointAndOffsets_test()
{

  auto const tol = boost::test_tools::tolerance(0.001);

  auto const [p, ofsA, ofsB] =
    geo::LineClosestPointAndOffsets(geo::origin() - 3.0 * geo::Xaxis(),
                                    geo::Xaxis(),
                                    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(),
                                    geo::Yaxis());

  BOOST_TEST(p.X() == 0.0, tol);
  BOOST_TEST(p.Y() == 0.0, tol);
  BOOST_TEST(p.Z() == 0.0, tol);
  BOOST_TEST(ofsA == +3.0, tol);
  BOOST_TEST(ofsB == -2.0, tol);

} // LineClosestPointAndOffsets_test()

// -----------------------------------------------------------------------------
void LineClosestPointWithScaledDirs_test()
{

  auto const tol = boost::test_tools::tolerance(0.001);

  auto const [p, ofsA, ofsB] =
    geo::LineClosestPointAndOffsets(geo::origin() - 3.0 * geo::Xaxis(),
                                    2.0 * geo::Xaxis(),
                                    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(),
                                    -2.0 * geo::Yaxis());

  BOOST_TEST(p.X() == 0.0, tol);
  BOOST_TEST(p.Y() == 0.0, tol);
  BOOST_TEST(p.Z() == 0.0, tol);
  BOOST_TEST(ofsA == (+3.0 / 2.0), tol);
  BOOST_TEST(ofsB == (-2.0 / -2.0), tol);

} // LineClosestPointWithScaledDirs_test()

// -----------------------------------------------------------------------------
void LineClosestPointWithNonHomogeneousDirs_test()
{

  auto const tol = boost::test_tools::tolerance(0.001);

  auto const [p, ofsA, ofsB] =
    geo::LineClosestPointAndOffsets(geo::origin() - 3.0 * geo::Xaxis(),
                                    1.5 * geo::Xaxis(),
                                    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(),
                                    -2.0 * geo::Yaxis());

  BOOST_TEST(p.X() == 0.0, tol);
  BOOST_TEST(p.Y() == 0.0, tol);
  BOOST_TEST(p.Z() == 0.0, tol);
  BOOST_TEST(ofsA == (+3.0 / 1.5), tol);
  BOOST_TEST(ofsB == (-2.0 / -2.0), tol);

} // LineClosestPointWithNonHomogeneousDirs_test()

// -----------------------------------------------------------------------------
void LineClosestPointAndOffsetsDocumentation_test()
{

  auto const tol = boost::test_tools::tolerance(0.001);

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
   * auto const xsectAndOfs = geo::LineClosestPointAndOffsets(
   *   geo::Point_t{ 0, 1, 0 }, geo::Vector_t{ 0.866, 0.0, 0.0 },
   *   geo::Point_t{ 2, 0, 1 }, geo::Vector_t{ 0.0,   0.5, 0.0 }
   *   );
   * point = xsectAndOfs.point;
   * offsetA = xsectAndOfs.offset1;
   * offsetB = xsectAndOfs.offset2;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (`point` to `geo::Point{ 2, 1, 0 }`, `offsetA` to `2.039...` and `offsetB`
   * to `2`, because the intersection point is always on the first line).
   * --- 8< --------------------------------------------------------------------
   */

  auto [point, offsetA, offsetB] = geo::LineClosestPointAndOffsets(geo::Point_t{2, 0, 1},
                                                                   geo::Vector_t{0.0, 0.5, 0.0},
                                                                   geo::Point_t{0, 1, 0},
                                                                   geo::Vector_t{0.866, 0.0, 0.0});

  // a way to check we did not mess with the assignment above too much
  static_assert(std::is_same_v<decltype(point), geo::Point_t>, "Unexpected point type");
  static_assert(std::is_same_v<decltype(offsetA), double>, "Unexpected first offset type");
  static_assert(std::is_same_v<decltype(offsetB), double>, "Unexpected second offset type");

  BOOST_TEST(point.X() == 2.0, tol);
  BOOST_TEST(point.Y() == 1.0, tol);
  BOOST_TEST(point.Z() == 1.0, tol);
  BOOST_TEST(offsetA == 2.0, tol);
  BOOST_TEST(offsetB == (2.0 / 0.866), tol);

  auto const xsectAndOfs = geo::LineClosestPointAndOffsets(geo::Point_t{0, 1, 0},
                                                           geo::Vector_t{0.866, 0.0, 0.0},
                                                           geo::Point_t{2, 0, 1},
                                                           geo::Vector_t{0.0, 0.5, 0.0});
  point = xsectAndOfs.point;
  offsetA = xsectAndOfs.offset1;
  offsetB = xsectAndOfs.offset2;

  BOOST_TEST(point.X() == 2.0, tol);
  BOOST_TEST(point.Y() == 1.0, tol);
  BOOST_TEST(point.Z() == 0.0, tol);
  BOOST_TEST(offsetA == 2.0 / 0.866, tol);
  BOOST_TEST(offsetB == 2.0, tol);

  // actually we _can_ assign with `std::tie()`:
  std::tie(point, offsetA, offsetB) =
    geo::LineClosestPointAndOffsets(geo::Point_t{2, 0, 1},
                                    geo::Vector_t{0.0, 0.5, 0.0},
                                    geo::Point_t{0, 1, 0},
                                    geo::Vector_t{0.866, 0.0, 0.0});

  BOOST_TEST(point.X() == 2.0, tol);
  BOOST_TEST(point.Y() == 1.0, tol);
  BOOST_TEST(point.Z() == 1.0, tol);
  BOOST_TEST(offsetA == 2.0, tol);
  BOOST_TEST(offsetB == (2.0 / 0.866), tol);

} // LineClosestPointAndOffsetsDocumentation_test()

// =============================================================================
void LineClosestPointWithUnitVectorsSimple_test()
{

  auto const tol = boost::test_tools::tolerance(0.001);

  geo::Point_t const p =
    geo::LineClosestPointWithUnitVectors(geo::origin() - geo::Zaxis() - 3.0 * geo::Xaxis(),
                                         geo::Xaxis(),
                                         geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(),
                                         geo::Yaxis());

  BOOST_TEST(p.X() == 0.0, tol);
  BOOST_TEST(p.Y() == 0.0, tol);
  BOOST_TEST(p.Z() == -1.0, tol);

} // LineClosestPointWithUnitVectorsSimple_test()

// -----------------------------------------------------------------------------
void LineClosestPointWithUnitVectorsSimple45_test()
{

  auto const tol = boost::test_tools::tolerance(0.001);

  geo::Point_t const p =
    geo::LineClosestPointWithUnitVectors(geo::origin(),
                                         geo::Xaxis(),
                                         geo::origin() + geo::Yaxis(),
                                         (geo::Xaxis() + geo::Zaxis()) / std::sqrt(2.0));

  BOOST_TEST(p.X() == 0.0, tol);
  BOOST_TEST(p.Y() == 0.0, tol);
  BOOST_TEST(p.Z() == 0.0, tol);

} // LineClosestPointWithUnitVectorsSimple45_test()

// -----------------------------------------------------------------------------
void LineClosestPointWithUnitVectorsAndOffsets_test()
{

  auto const tol = boost::test_tools::tolerance(0.001);

  auto const [p, ofsA, ofsB] = geo::LineClosestPointAndOffsetsWithUnitVectors(
    geo::origin() - 3.0 * geo::Xaxis(),
    geo::Xaxis(),
    geo::origin() + geo::Zaxis() + 2.0 * geo::Yaxis(),
    geo::Yaxis());

  BOOST_TEST(p.X() == 0.0, tol);
  BOOST_TEST(p.Y() == 0.0, tol);
  BOOST_TEST(p.Z() == 0.0, tol);
  BOOST_TEST(ofsA == +3.0, tol);
  BOOST_TEST(ofsB == -2.0, tol);

} // LineClosestPointWithUnitVectorsAndOffsets_test()

// -----------------------------------------------------------------------------
void LineClosestPointAndOffsetsWithUnitVectorsDocumentation_test()
{

  auto const tol = boost::test_tools::tolerance(0.001);

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
   * auto const xsectAndOfs = geo::LineClosestPointAndOffsetsWithUnitVectors(
   *   geo::Point_t{ 0, 1, 0 }, geo::Vector_t{ 1, 0, 0 },
   *   geo::Point_t{ 2, 0, 1 }, geo::Vector_t{ 0, 1, 0 }
   *   );
   * point = xsectAndOfs.point;
   * offsetA = xsectAndOfs.offset1;
   * offsetB = xsectAndOfs.offset2;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (`point` to `geo::Point{ 2, 1, 0 }`, `offsetA` to `2` and `offsetB` to `1`,
   * because the intersection point is always on the first line).
   * --- 8< --------------------------------------------------------------------
   */

  auto [point, offsetA, offsetB] = geo::LineClosestPointAndOffsetsWithUnitVectors(
    geo::Point_t{2, 0, 1}, geo::Vector_t{0, 1, 0}, geo::Point_t{0, 1, 0}, geo::Vector_t{1, 0, 0});

  // a way to check we did not mess with the assignment above too much
  static_assert(std::is_same_v<decltype(point), geo::Point_t>, "Unexpected point type");
  static_assert(std::is_same_v<decltype(offsetA), double>, "Unexpected first offset type");
  static_assert(std::is_same_v<decltype(offsetB), double>, "Unexpected second offset type");

  BOOST_TEST(point.X() == 2.0, tol);
  BOOST_TEST(point.Y() == 1.0, tol);
  BOOST_TEST(point.Z() == 1.0, tol);
  BOOST_TEST(offsetA == 1.0, tol);
  BOOST_TEST(offsetB == 2.0, tol);

  auto const xsectAndOfs = geo::LineClosestPointAndOffsetsWithUnitVectors(
    geo::Point_t{0, 1, 0}, geo::Vector_t{1, 0, 0}, geo::Point_t{2, 0, 1}, geo::Vector_t{0, 1, 0});
  point = xsectAndOfs.point;
  offsetA = xsectAndOfs.offset1;
  offsetB = xsectAndOfs.offset2;

  BOOST_TEST(point.X() == 2.0, tol);
  BOOST_TEST(point.Y() == 1.0, tol);
  BOOST_TEST(point.Z() == 0.0, tol);
  BOOST_TEST(offsetA == 2.0, tol);
  BOOST_TEST(offsetB == 1.0, tol);

  // actually we _can_ assign with `std::tie()`:
  std::tie(point, offsetA, offsetB) = geo::LineClosestPointAndOffsetsWithUnitVectors(
    geo::Point_t{2, 0, 1}, geo::Vector_t{0, 1, 0}, geo::Point_t{0, 1, 0}, geo::Vector_t{1, 0, 0});

  BOOST_TEST(point.X() == 2.0, tol);
  BOOST_TEST(point.Y() == 1.0, tol);
  BOOST_TEST(point.Z() == 1.0, tol);
  BOOST_TEST(offsetA == 1.0, tol);
  BOOST_TEST(offsetB == 2.0, tol);

} // LineClosestPointAndOffsetsWithUnitVectorsDocumentation_test()

// =============================================================================
BOOST_AUTO_TEST_CASE(LineClosestPointTestCase)
{

  LineClosestPointSimple_test();
  LineClosestPointSimple45_test();
  LineClosestPointAndOffsets_test();
  LineClosestPointWithScaledDirs_test();
  LineClosestPointWithNonHomogeneousDirs_test();
  LineClosestPointAndOffsetsDocumentation_test();

} // BOOST_AUTO_TEST_CASE(LineClosestPointTestCase)

// -----------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(LineClosestPointWithUnitVectorsTestCase)
{

  LineClosestPointWithUnitVectorsSimple_test();
  LineClosestPointWithUnitVectorsSimple45_test();
  LineClosestPointWithUnitVectorsAndOffsets_test();
  LineClosestPointAndOffsetsWithUnitVectorsDocumentation_test();

} // BOOST_AUTO_TEST_CASE(LineClosestPointWithUnitVectorsTestCase)

// -----------------------------------------------------------------------------
