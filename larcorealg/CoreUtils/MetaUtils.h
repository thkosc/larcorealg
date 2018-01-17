/**
 * @addtogroup Utilities General utilities
 * @brief General programming utilities.
 * 
 * @{
 */
/**
 * @defgroup Metaprogramming General utilities for metaprogramming
 * @brief General utilities for use with templates and metaprogramming.
 */
/**
 * @}
 */

/** ****************************************************************************
 * @file   larcorealg/CoreUtils/MetaUtils.h
 * @brief  Basic C++ metaprogramming utilities.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   January 12, 2017
 * @ingroup Metaprogramming
 */

#ifndef LARCOREALG_COREUTILS_METAUTILS_H
#define LARCOREALG_COREUTILS_METAUTILS_H


// C/C++ standard libraries
#include <type_traits>


/**
 * @namespace util
 * @brief Namespace for general, non-LArSoft-specific utilities.
 * @ingroup Utilities
 */
namespace util {
  
  //--- BEGIN MetaprogrammingBase ----------------------------------------------
  /**
   * @defgroup MetaprogrammingBase Simple utility traits
   * @brief Simple traits for the implementation of other traits.
   * @ingroup Metaprogramming
   */
  /// @{
  
  //----------------------------------------------------------------------------
  /// Trait returning the very same type as in the template argument.
  template <typename T>
  struct self_type { using type = T; };
  
  /// The very same type as in the template argument.
  template <typename T>
  using self_t = typename self_type<T>::type;
  
  
  //----------------------------------------------------------------------------
  /**
   * @brief A `std::false_type` with a template argument.
   * @see util::always_true_type, util::always_false_v
   * 
   * This type allows a `static_assert` to fail only when the template type it's
   * in is being instantiated:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T>
   * struct MandatoryCustomizationPoint {
   *   static_assert(util::always_false_type<T>(), "You have to specialize!");
   * };
   * 
   * template <typename T>
   * struct MandatoryCustomizationPoint<std::vector<T>> {
   *   using type = typename std::vector<T>::reference;
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In this example, using `std::false_type` might have tipped the compiler to
   * trigger the assertion failure even if the class is not instantiated.
   */
  template <typename>
  struct always_false_type: public std::false_type {};
  
  
  //----------------------------------------------------------------------------
  /**
   * @brief A templated constant, always false.
   * @see util::always_false_type, util::always_true_v
   * 
   * This constant allows a `static_assert` to fail only when the template type
   * it's in is being instantiated:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T>
   * struct MandatoryCustomizationPoint {
   *   static_assert(util::always_false_v<T>, "You have to specialize!");
   * };
   * 
   * template <typename T>
   * struct MandatoryCustomizationPoint<std::vector<T>> {
   *   using type = typename std::vector<T>::reference;
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In this example, using `std::false_type` might have tipped the compiler to
   * trigger the assertion failure even if the class is not instantiated.
   */
  template <typename>
  constexpr bool always_false_v = false;
  
  
  /**
   * @brief A `std::true_type` with a template argument.
   * @see util::always_false_type, util::always_true_v
   * 
   * This is one way to allow to specialize for classes with a certain type:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T, typename = void>
   * class ReferenceTypeExtractor {
   *   static_assert(util::always_false_type<T>(), "Type has no reference!");
   * };
   * 
   * template <typename Cont>
   * struct ReferenceTypeExtractor<
   *   Cont,
   *   std::enable_if_t
   *     <util::always_true_type<typename Cont::value_type>::value>
   *   >
   * {
   *   using type = typename Cont::reference;
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename>
  struct always_true_type: public std::true_type {};
  
  
  /**
   * @brief A template constant always true.
   * @see util::always_true_type, util::always_false_v
   * 
   * This is one way to allow to specialize for classes with a certain type:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * template <typename T, typename = void>
   * class ReferenceTypeExtractor {
   *   static_assert(util::always_false_v<T>, "Type has no reference!");
   * };
   * 
   * template <typename Cont>
   * struct ReferenceTypeExtractor
   *   <Cont, std::enable_if_t<util::always_true_v<typename Cont::value_type>>>
   * {
   *   using type = typename Cont::reference;
   * };
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   */
  template <typename>
  constexpr bool always_true_v = true;
  
  
  /// @}
  //--- END MetaprogrammingBase ------------------------------------------------
  
  
} // namespace util


#endif // LARCOREALG_COREUTILS_METAUTILS_H

