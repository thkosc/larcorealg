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
  
  namespace details {
    
    /// Implementation detail of `staticDumpClassName()`.
    template <typename T>
    struct ClassNameStaticDumper;
    
  } // namespace details
  
  
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
  
  
  // will be present in C++17 as std::bool_constant
  /// @todo Replace all uses with std::bool_constant when C++17 is adopted
  template <bool Value>
  using bool_constant = std::integral_constant<bool, Value>;
  
  // will be present in C++17 as std::negation
  /// @todo Replace all uses with std::negation when C++17 is adopted
  template <typename BoolTrait>
  using negation = bool_constant<!BoolTrait::value>;
  
  /// The negation of `std::is_same`.
  /// @todo Switch to std::negation in implementation when C++17 is adopted
  template <typename A, typename B>
  using is_not_same = negation<std::is_same<A, B>>;
  
  
  // @{
  /**
   * @brief Helper to determine the type of a variable at compilation time.
   * @tparam T type to be investigated
   * 
   * It may be difficult to understand which type is being used in a failing
   * static assertion or in some complicate metaprogramming code (is there any
   * other kind?), especially when due to a compilation failure the code can't
   * be run.
   * This class is supposed to help by forcing the compiler to a halt, and it is
   * devised so that the compiler should print in the error message what it
   * thinks the type `T` is.
   * 
   * An example of usage:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * #include "MetaUtils.h"
   * 
   * void f() {
   *   constexpr auto v = 5U - 6U; // which type is `v` of?
   *   util::staticDumpClassName(v);
   * }
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * For example, Clang 5.0.1 emits these errors:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In file included from metatest.cpp:1:
   * ./MetaUtils.h:217:7: error: static_assert failed "ClassNameStaticDumper<T>: look for T in the error message context"
   *       static_assert(
   *       ^
   * ./MetaUtils.h:228:39: note: in instantiation of template class 'util::details::ClassNameStaticDumper<unsigned int>' requested here
   *   void staticDumpClassName() { (void) details::ClassNameStaticDumper<T>(); }
   *                                       ^
   * ./MetaUtils.h:195:46: note: in instantiation of function template specialization 'util::staticDumpClassName<unsigned int>' requested here
   *   [[noreturn]] void staticDumpClassName(T) { staticDumpClassName<T>(); }
   *                                              ^
   * metatest.cpp:5:16: note: in instantiation of function template specialization 'util::staticDumpClassName<unsigned int>' requested here
   *   (void) util::staticDumpClassName(v);
   *                ^
   * 1 error generated.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * From the first "note" we can see that the type of `v` is `unsigned int`.
   * Note that if `v` is not copiable, an additional error will be emitted.
   * To avoid that, the type can be specified as template parameter, as in
   * `util::staticDumpClassName<decltype(v)>()`.
   * The same program is processed by GNU GCC with an error message:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * In file included from metatest.cpp:1:0:
   * MetaUtils.h: In instantiation of ‘struct util::details::ClassNameStaticDumper<unsigned int>’:
   * MetaUtils.h:248:48:   required from ‘void util::staticDumpClassName() [with T = unsigned int]’
   * MetaUtils.h:215:68:   required from ‘void util::staticDumpClassName(T) [with T = unsigned int]’
   * metatest.cpp:5:30:   required from here
   * MetaUtils.h:237:7: error: static assertion failed: ClassNameStaticDumper<T>: look for T in the error message context
   *        static_assert(
   *        ^~~~~~~~~~~~~
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * where the type is mentioned in all the three context lines.
   */
  template <typename T>
  void staticDumpClassName();
  
  template <typename T>
  void staticDumpClassName(T) { staticDumpClassName<T>(); }
  // @}
  
  
  /// @}
  //--- END MetaprogrammingBase ------------------------------------------------
  
  
} // namespace util


//------------------------------------------------------------------------------
//---  Template implementation
//------------------------------------------------------------------------------
namespace util {
  
  //----------------------------------------------------------------------------
  namespace details {
    
    /// Implementation detail of `staticDumpClassName()`.
    template <typename T>
    struct ClassNameStaticDumper {
      static_assert(
        util::always_false_type<T>(),
        "ClassNameStaticDumper<T>: look for T in the error message context"
        );
    }; // struct ClassNameStaticDumper
    
  } // namespace details
  
  
  //----------------------------------------------------------------------------
  template <typename T>
  void staticDumpClassName() { (void) details::ClassNameStaticDumper<T>(); }
  
  
  //----------------------------------------------------------------------------
  
} // namespace util


//------------------------------------------------------------------------------

#endif // LARCOREALG_COREUTILS_METAUTILS_H

