#include "larcorealg/Geometry/Intersections.h"
#include "larcorealg/CoreUtils/NumericUtils.h"
#include "larcorealg/CoreUtils/RealComparisons.h"

#include "messagefacility/MessageLogger/MessageLogger.h"

namespace geo {

  // Functions to allow determination if two wires intersect, and if so where.
  // This is useful information during 3D reconstruction.
  //......................................................................
  bool IntersectLines(double A_start_x,
                      double A_start_y,
                      double A_end_x,
                      double A_end_y,
                      double B_start_x,
                      double B_start_y,
                      double B_end_x,
                      double B_end_y,
                      double& x,
                      double& y)
  {
    // Equation from http://en.wikipedia.org/wiki/Line%E2%80%93line_intersection
    // T.Yang Nov, 2014
    // Notation: x => coordinate orthogonal to the drift direction and to the beam direction
    //           y => direction orthogonal to the previous and to beam direction

    double const denom =
      (A_start_x - A_end_x) * (B_start_y - B_end_y) - (A_start_y - A_end_y) * (B_start_x - B_end_x);

    constexpr lar::util::RealComparisons<double> coordIs{1e-8};
    if (coordIs.zero(denom)) return false;

    double const A = (A_start_x * A_end_y - A_start_y * A_end_x) / denom;
    double const B = (B_start_x * B_end_y - B_start_y * B_end_x) / denom;

    x = (B_start_x - B_end_x) * A - (A_start_x - A_end_x) * B;
    y = (B_start_y - B_end_y) * A - (A_start_y - A_end_y) * B;

    return true;
  }

  //......................................................................
  bool IntersectSegments(double A_start_x,
                         double A_start_y,
                         double A_end_x,
                         double A_end_y,
                         double B_start_x,
                         double B_start_y,
                         double B_end_x,
                         double B_end_y,
                         double& x,
                         double& y)
  {
    bool const bCross = IntersectLines(
      A_start_x, A_start_y, A_end_x, A_end_y, B_start_x, B_start_y, B_end_x, B_end_y, x, y);

    if (bCross) {
      mf::LogWarning("IntersectSegments") << "The segments are parallel!";
      return false;
    }

    return lar::util::PointWithinSegments(
      A_start_x, A_start_y, A_end_x, A_end_y, B_start_x, B_start_y, B_end_x, B_end_y, x, y);
  }

} // namespace geo
