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
#include <array>
#include <functional> // std::reference_wrapper<>
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
  
  
  //--- BEGIN Type identification ----------------------------------------------
  /**
   * @defgroup MetaprogrammingTypeIdentification Determination of specific types
   * @brief Traits to identify specific types.
   * @ingroup Metaprogramming
   */
  /// @{
  
  //----------------------------------------------------------------------------
  /**
   * @brief Identifies whether the specified type is a STL array.
   * @tparam T the type to be tested
   * @see `util::is_STLarray_v`
   */
  template <typename>
  struct is_STLarray: public std::false_type {};
  
  /**
   * @brief A constant describing whether the specified type is a STL array.
   * @tparam T the type to be tested
   * @see `util::is_STLarray`
   */
  template <typename T>
  constexpr bool is_STLarray_v = is_STLarray<T>::value;
  
  
  //----------------------------------------------------------------------------
  /**
   * @brief Identifies whether the specified type is a `std::reference_wrapper`.
   * @tparam T the type to be tested
   * @see `util::is_reference_wrapper_v`
   */
  template <typename T>
  struct is_reference_wrapper;
  
  /**
   * @brief A constant describing whether the specified type is a
   *        `std::reference_wrapper`.
   * @tparam T the type to be tested
   * @see `util::is_STLarray`
   */
  template <typename T>
  constexpr bool is_reference_wrapper_v = is_reference_wrapper<T>::value;
  
  
  /// @}
  //--- END Type identification ------------------------------------------------
  
  
  //--- BEGIN Type manipulation ------------------------------------------------
  /**
   * @defgroup MetaprogrammingTypeManipulation Manipulation of types
   * @brief Traits to change types.
   * @ingroup Metaprogramming
   */
  /// @{
  
  //----------------------------------------------------------------------------
  /**
   * @brief Trait with type `Base`, plus the constantness as in `Key`.
   * @tparam Base the basic type being returned
   * @tparam Key a type expressing the constantness wanted for `Base`
   * 
   * The `type` member of this trait is:
   * * the type `Base` with constantness removed, if `Key` is non-constant type
   * * the type `Base` with constantness added, if `Key` is a constant type
   * 
   * This trait passes through references. Both `Base` and `Key` are treated as
   * they were no references (e.g. a `int const&` is treated as a `int const`).
   * The referenceness of the type in the trait is the same as the one of
   * `Base`.
   * 
   * Therefore, for example:
   * * `with_const_as<int const, double>` yields `int`;
   * * `with_const_as<int const, double const>` yields `int const`;
   * * `with_const_as<int const&, double>` yields `int&`;
   * * `with_const_as<int, double const&>` yields `int const`.
   * 
   * R-value references are likewise taken into account.
   */
  template <typename Base, typename Key>
  struct with_const_as;
  
  /**
   * @brief The type `Base`, plus the constantness as in `Key`.
   * @tparam Base the basic type being returned
   * @tparam Key a type expressing the constantness wanted for `Base`
   * @see `util::with_const_as`
   */
  template <typename Base, typename Key>
  using with_const_as_t = typename with_const_as<Base, Key>::type;
  
  
  //----------------------------------------------------------------------------
  /**
   * @brief Trait with type `T` stripped of all known reference types.
   * @tparam T type to remove referenceness from
   * 
   * In addition of the standard C++ references, this trait also removes all
   * pseudo-referenceness known (to it).
   * That currently includes:
   * * C++ l-value and r-value reference
   * * `std::reference_wrapper`
   * 
   * The implementation removes all references recursively.
   */
  template <typename T>
  struct strip_referenceness_type;
  
  /**
   * @brief The type `T` stripped of all known reference types.
   * @tparam T type to remove referenceness from
   * @see `util::strip_referenceness_type`
   */
  template <typename T>
  using strip_referenceness_t = typename strip_referenceness_type<T>::type;
  
  
  //----------------------------------------------------------------------------
  
  
  /// @}
  //--- END Type manipulation --------------------------------------------------
  
} // namespace util


//------------------------------------------------------------------------------
//---  Template implementation
//------------------------------------------------------------------------------
namespace util {
  
  //----------------------------------------------------------------------------
  namespace details {
    
    //--------------------------------------------------------------------------
    /// Implementation detail of `staticDumpClassName()`.
    template <typename T>
    struct ClassNameStaticDumper {
      static_assert(
        util::always_false_type<T>(),
        "ClassNameStaticDumper<T>: look for T in the error message context"
        );
    }; // struct ClassNameStaticDumper
    
    //--------------------------------------------------------------------------
    // implementation for `with_const_as`
    
    // - final implementation:
    template <typename Base, typename /* Key */, typename = void>
    struct with_const_as_impl
      { using type = std::remove_const_t<Base>; };
    
    template <typename Base, typename Key>
    struct with_const_as_impl
      <Base, Key, std::enable_if_t<std::is_const_v<Key>>>
      { using type = std::add_const_t<Base>; };
    
    // - implementation dispatcher for reference types
    //   - pass through for not-reference types
    template <typename Base, typename Key, typename = void>
    struct with_const_as_dispatch_ref: with_const_as_impl<Base, Key> {};
    //   - lvalue reference
    template <typename Base, typename Key>
    struct with_const_as_dispatch_ref
      <Base, Key, std::enable_if_t<std::is_lvalue_reference_v<Base>>>
    {
      using type = std::add_lvalue_reference_t
        <typename with_const_as_impl<std::remove_reference_t<Base>, Key>::type>;
    };
    //   - rvalue reference
    template <typename Base, typename Key>
    struct with_const_as_dispatch_ref
      <Base, Key, std::enable_if_t<std::is_rvalue_reference_v<Base>>>
    {
      using type = std::add_rvalue_reference_t
        <typename with_const_as_impl<std::remove_reference_t<Base>, Key>::type>;
    };
    
    // - key management
    template <typename Base, typename Key>
    struct with_const_as_dispatch_keyref
      : with_const_as_dispatch_ref<Base, std::remove_reference_t<Key>> {};
    
    // - top level implementation dispatcher
    template <typename Base, typename Key>
    struct with_const_as_dispatcher: with_const_as_dispatch_keyref<Base, Key> {};
    
    
    //--------------------------------------------------------------------------
    //--- implementation of `is_reference_wrapper`
    template <typename T>
    struct is_reference_wrapper_impl: std::false_type {};
    
    template <typename T>
    struct is_reference_wrapper_impl<std::reference_wrapper<T>>
      : public std::true_type
    {};
  
  
    //--------------------------------------------------------------------------
    //--- implementation of `strip_referenceness_type`
    
    template <typename T>
    struct strip_referenceness_type_impl;
    
    // implementation layer dealing with `std::reference_wrapper`
    template <typename T, typename = void>
    struct strip_referenceness_type_impl_wrapref
      { using type = T; }; // exit here
    
    // - handle any constantness and volatility
    template <typename T>
    struct strip_referenceness_type_impl_wrapref<
      T,
      std::enable_if_t<util::is_reference_wrapper_v<std::remove_cv_t<T>>>
      >
      : strip_referenceness_type_impl<typename T::type> // back to square one
    {}; 
    
    // implementation layer dealing with C++ references
    template <typename T>
    struct strip_referenceness_type_impl_ref
      : strip_referenceness_type_impl_wrapref<T> {};
    
    template <typename T>
    struct strip_referenceness_type_impl_ref<T&>
      : strip_referenceness_type_impl_wrapref<T> {};
    
    template <typename T>
    struct strip_referenceness_type_impl_ref<T&&>
      : strip_referenceness_type_impl_wrapref<T> {};
    
    // entry point: start by dealing with C++ references
    template <typename T>
    struct strip_referenceness_type_impl: strip_referenceness_type_impl_ref<T>
    {};
    
    
    //--------------------------------------------------------------------------
    
    
  } // namespace details
  
  
  //----------------------------------------------------------------------------
  template <typename T>
  void staticDumpClassName() { (void) details::ClassNameStaticDumper<T>(); }
  
  
  //----------------------------------------------------------------------------
  template <typename T, std::size_t N>
  struct is_STLarray<std::array<T, N>>: public std::true_type {};
  
  
  //----------------------------------------------------------------------------
  template <typename T>
  struct is_reference_wrapper:
    details::is_reference_wrapper_impl<std::decay_t<std::remove_reference_t<T>>>
    {};
  
  
  //----------------------------------------------------------------------------
  template <typename Base, typename Key>
  struct with_const_as: public details::with_const_as_dispatcher<Base, Key> {};
  
  //----------------------------------------------------------------------------
  template <typename T>
  struct strip_referenceness_type
    : public details::strip_referenceness_type_impl<T>
  {};
  
  //----------------------------------------------------------------------------
  
} // namespace util


//------------------------------------------------------------------------------

#endif // LARCOREALG_COREUTILS_METAUTILS_H

