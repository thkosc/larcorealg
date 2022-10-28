#ifndef LARCOREALG_GEOMETRY_INTERSECTIONS_H
#define LARCOREALG_GEOMETRY_INTERSECTIONS_H

namespace geo {
  // Functions to allow determination if two wires intersect, and if so where.
  // This is useful information during 3D reconstruction.
  //......................................................................

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

#endif // LARCOREALG_GEOMETRY_INTERSECTIONS_H
