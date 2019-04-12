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
#include <type_traits>


namespace util {

  // @{
  /**
   * @brief Returns the absolute value of the difference between two values.
   * @tparam A type of the first value
   * @tparam B type of the second value (*must* actually be as `A`)
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
   * It still assumes that the difference is representable in `A`; for example,
   * this assumption will fail for `int` types with `a` a very large number and
   * `b` a very small (i.e. negative) number.
   *
   * Requirements:
   * * `A` and `B` must be the same type
   *
   */
  template <typename A, typename B>
  constexpr auto absDiff(A const& a, B const& b)
    {
      static_assert(
        std::is_same<std::decay_t<A>, std::decay_t<B>>(),
        "Arguments of util::absDiff() have to be of the same type."
        );
      return (b > a)? (b - a): (a - b);
    }
  // @}

} // namespace util


#endif // LARCOREALG_COREUTILS_NUMERICUTILS_H

