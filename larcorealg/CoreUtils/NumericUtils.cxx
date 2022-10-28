#include "larcorealg/CoreUtils/NumericUtils.h"
#include "larcorealg/CoreUtils/RealComparisons.h"

#include <cmath>

namespace lar::util {
  bool ValueInRange(double value, double min, double max)
  {
    if (min > max) std::swap(min, max); //protect against funny business due to wire angles
    if (std::abs(value - min) < 1e-6 || std::abs(value - max) < 1e-6) return true;
    return (value >= min) && (value <= max);
  }

  //--------------------------------------------------------------------
  /// Returns whether x and y are within both specified ranges (A and B).
  bool PointWithinSegments(double A_start_x,
                           double A_start_y,
                           double A_end_x,
                           double A_end_y,
                           double B_start_x,
                           double B_start_y,
                           double B_end_x,
                           double B_end_y,
                           double x,
                           double y)
  {
    constexpr RealComparisons<double> coordIs{1e-8};
    return coordIs.withinSorted(x, A_start_x, A_end_x) &&
           coordIs.withinSorted(y, A_start_y, A_end_y) &&
           coordIs.withinSorted(x, B_start_x, B_end_x) &&
           coordIs.withinSorted(y, B_start_y, B_end_y);
  }

}
