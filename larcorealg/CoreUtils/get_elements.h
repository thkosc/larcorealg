/**
 * @file   larcorealg/CoreUtils/get_elements.h
 * @brief  Definition of `util::get_elements()` and `util::get_const_elements()`.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   December 13, 2019
 *
 * This is a header-only library.
 */

#ifndef LARCOREALG_COREUTILS_GET_ELEMENTS_H
#define LARCOREALG_COREUTILS_GET_ELEMENTS_H

// LArSoft libraries
#include "larcorealg/CoreUtils/StdUtils.h" // util::get()
#include "larcorealg/CoreUtils/span.h"     // util::make_transformed_span()

// C/C++ libraries
#include <cstddef> // std::size_t
#include <type_traits>
#include <utility> // std::get()

namespace util {

  // -- BEGIN -- Transformed iterations ----------------------------------------
  /// @name Transformed iterations
  /// @{

  /**
   * @brief Range-for loop helper iterating across some of the element of
   *        each value in the specified collection.
   * @tparam Indices indices of the elements to extract
   * @tparam Coll type of the collection to iterate through
   * @param coll the collection to iterate through
   * @return an object suitable for range-for loop
   * @see `util::get_const_elements()`
   *
   * This function enables range-for loops with a selection of elements from
   * a container of `tuple` (or anything responding to `util::get()`).
   *
   * The following example fills a map using as key the first element (`0U`)
   * of a tuple and as value the third element (`2U`):
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<std::tuple<char, int, float>> data
   *   { { 'z', 0, 1.0F }, { 'o', 1, 3.0F }, { 't', 2, 9.0F } };
   * std::map<char, double> factors;
   *
   * for (auto const& [ letter, factor ]: util::get_elements<0U, 2U>(data)) {
   *   factors.emplace(letter, factor);
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * If only one index is specified, the loop will not use structured binding,
   * but rather a simple variable:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * std::vector<int> exponents;
   *
   * for (int exponent: util::get_elements<1U>(data)) {
   *   exponents.push_back(exponent);
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *
   * While the examples do not demonstrate changing the values in the `data`
   * collection, that is also supported.
   *
   * @note This function also works with associative containers based on
   *       `std::pair` (`std::map` and the likes).
   */
  template <std::size_t... Indices, typename Coll>
  decltype(auto) get_elements(Coll&& coll);

  /**
   * @brief Range-for loop helper iterating across the constant values of the
   *        specified collection.
   * @see `util::values()`
   *
   * This function is equivalent to `util::values()` but the values are
   * extracted as if the specified collection were constant.
   */
  template <std::size_t... Indices, typename Coll>
  decltype(auto) get_const_elements(Coll&& coll);

  /// @}
  // -- END -- Transformed iterations ------------------------------------------

} // namespace util

//==============================================================================
//=== template implementation
//==============================================================================
//------------------------------------------------------------------------------
//--- util::values()
//------------------------------------------------------------------------------
namespace util::details {

  //----------------------------------------------------------------------------
  template <typename Coll, std::size_t First, std::size_t... Others>
  struct get_elements_impl {

    template <typename T>
    static constexpr decltype(auto) iterate(T&& coll) noexcept
    {
      auto extractor = [](auto&& item) -> decltype(auto) {
        if constexpr (sizeof...(Others) == 0U) { return util::get<First>(item); }
        else {
          return std::forward_as_tuple(util::get<First>(item), util::get<Others...>(item));
        }
      };
      return util::make_transformed_span(coll, extractor);
    }

  }; // struct get_elements_impl

  //----------------------------------------------------------------------------

} // namespace util::details

//------------------------------------------------------------------------------
template <std::size_t... Indices, typename Coll>
decltype(auto) util::get_elements(Coll&& coll)
{
  static_assert(sizeof...(Indices) >= 1U, "At least one index must be specified");
  return details::get_elements_impl<std::decay_t<Coll>, Indices...>::iterate(
    std::forward<Coll>(coll));
} // util::get_elements()

//------------------------------------------------------------------------------
template <std::size_t... Indices, typename Coll>
decltype(auto) util::get_const_elements(Coll&& coll)
{
  return get_elements<Indices...>(std::as_const(coll));
}

//------------------------------------------------------------------------------

#endif // LARCOREALG_COREUTILS_GET_ELEMENTS_H
