/**
 * @file   larcorealg/CoreUtils/NumericUtils.h
 * @brief  Functions to help with numbers.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   March 30, 2018
 * 
 * This library is currently header-only.
 */
#ifndef LARCOREALG_COREUTILS_NUMERICUTILS_H
#define LARCOREALG_COREUTILS_NUMERICUTILS_H

// C/C++ standard libraries

namespace util {
  
  // @{
  /**
   * @brief Returns the absolute value of the difference between two values.
   * @tparam A type of the first value
   * @tparam B type of the second value
   * @param a the first value
   * @param b the second value
   * @return the difference between the largest and the smallest of `a` and `b`
   * 
   * The pecularity of this implementation is that it always avoids taking the
   * difference between the smallest and the largest of `a` and `b`. An
   * equivalent implementation is:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * return std::max(a, b) - std::min(a, b);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * 
   * 
   * Requirements:
   * * `A` and `B` must be comparable: in particular, `operator< (B, A)` must
   *     exist
   * * it must be possible to subtract `A` and `B`: in particular, both
   *     `operator- (A, B)` and `operator- (B, A)` must exist and yield a result
   *     of the same data type
   * 
   * 
   */
  template <typename A, typename B>
  constexpr auto absDiff(A const& a, B const& b)
    { return (b > a)? (b - a): (a - b); }
  // @}
  
} // namespace util


#endif // LARCOREALG_COREUTILS_NUMERICUTILS_H

