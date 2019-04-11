/** ****************************************************************************
 * @file   larcorealg/CoreUtils/ContainerMeta.h
 * @brief  C++ metaprogramming utilities for dealing with containers.
 * @author Gianluca Petrillo (petrillo@fnal.gov)
 * @date   July 27, 2017
 * @ingroup Metaprogramming
 */

#ifndef LARCOREALG_COREUTILS_CONTAINERMETA_H
#define LARCOREALG_COREUTILS_CONTAINERMETA_H

// LArSoft 
#include "larcorealg/CoreUtils/MetaUtils.h" // for Doxygen documentation

// C/C++ standard libraries
#include <iterator> // std::begin(), std::cbegin()
#include <functional> // std::reference_wrapper<>
#include <memory> // std::unique_ptr<>
#include <utility> // std::declval()
#include <type_traits> // std::enable_if_t<>


namespace util {
  
  
  //--- BEGIN ContainerMetaprogramming -----------------------------------------
  /**
   * @defgroup ContainerMetaprogramming Traits for containers
   * @brief Simple traits for container classes.
   * @ingroup Metaprogramming
   * 
   * Trait classes describing a type have a `type` member defined after that
   * type. They are also available as template types, in the fashion of C++14.
   * 
   * Trait classes describing a value have a `value` member (static, constexpr)
   * defined with that value. They are also available as template constants,
   * in the fashion of C++17.
   * 
   * 
   * Container value traits
   * -----------------------
   * 
   * Container classes mimicking the C++ STL interface provide a `value_type`
   * type detailing the contained object. This is generalised in the
   * `util::collection_value_type` trait, which can be specialised for
   * containers which do not support that.
   * 
   * Two additional traits are provided, which describe the type obtained by
   * accessing an element of the container. This is a close relative of
   * `util::collection_value_type`, but it may be decorated by l-value reference
   * or constantness depending on the container. The two traits represent access
   * one to a container that is mutable, the other to one that is constant.
   * 
   */
  /// @{
  
  
  //----------------------------------------------------------------------------
  /// Trait of value contained in the template collection `Coll`.
  template <typename Coll>
  struct collection_value_type;
  
  /// Type contained in the collection `Coll`.
  template <typename Coll>
  using collection_value_t = typename collection_value_type<Coll>::type;
  
  
  //----------------------------------------------------------------------------
  /// Trait of type obtained by access to element of collection `Coll`.
  template <typename Coll>
  struct collection_value_access_type;
  
  /// Type obtained by constant access to element of collection `Coll`.
  template <typename Coll>
  using collection_value_access_t
    = typename collection_value_access_type<Coll>::type;
  
  
  //----------------------------------------------------------------------------
  /// Trait of type obtained by constant access to element of collection `Coll`.
  template <typename Coll>
  struct collection_value_constant_access_type;
  
  /// Type obtained by constant access to element of collection `Coll`.
  template <typename Coll>
  using collection_value_constant_access_t
    = typename collection_value_constant_access_type<Coll>::type;
  
  
  //----------------------------------------------------------------------------
  /**
   * @brief Trait of a type that can be used to reference the collection `Coll`.
   * @tparam Coll type of the collection to be referenced
   * 
   * The goal is to have an object with access to the data of the collection of
   * type `Coll`: this object should be able to be copied, but should not
   * duplicate (copy) the data.
   * The most versatile solution is to have a reference to `Coll`, and in
   * particular `std::reference_wrapper` does the job. But in some special cases
   * `Coll` itself will do already, as it is for a bare C pointer.
   * 
   * This trait describes a type with these characteristics, privileging the
   * simplest solution.
   */
  template <typename Coll>
  struct collection_reference_type;
  
  /// The type contained in `util::collection_reference_type` trait.
  template <typename Coll>
  using collection_reference_t = typename collection_reference_type<Coll>::type;
  
  
  /**
   * @brief Returns an object referencing to the data contained in `coll`.
   * @tparam Coll type of collection of the data
   * @return an object referencing to the data contained in `coll`
   * @see `util::collection_from_reference()`
   * 
   * The criteria, as well as the type of the returned object, are similar to
   * `util::collection_reference_type`.
   * Therefore, for example a C pointer is returned unchanged, while a
   * `std::vector` is returned wrapped into a `std::reference_wrapper`.
   * A `std::unique_ptr` is returned as bare pointer, given that the returned
   * object does not own the data.
   * 
   */
  template <typename Coll>
  auto make_collection_reference(Coll&& coll);
  
  
  /**
   * @brief Trait with the type of collection referenced by `collRef`.
   * @tparam CollRef type of collection of the data
   * @see `util::make_collection_reference()`,
   *      `util::collection_from_reference()`
   * 
   * The type is a direct reference to the unwrapper container.
   * For example, a `CollRef` instance of `std::reference_wrapper` will result
   * in a C reference of the wrapped container, while a C pointer is left
   * unchanged and a `std::unique_ptr` is turned into the equivalent pointer
   * to its elements.
   */
  template <typename Cont>
  struct collection_from_reference_type;
  
  /// Type contained in `util::collection_from_reference_type` trait.
  template <typename Cont>
  using collection_from_reference_t
    = typename collection_from_reference_type<Cont>::type;
  
  /**
   * @brief Returns the object referenced by `collRef` as a C++ reference.
   * @tparam CollRef type of collection of the data
   * @param collRef collection of the data to be referenced to
   * @return a reference to the object referenced by `collRef`
   * @see `util::make_collection_reference()`
   * 
   * The return value is a direct reference to the unwrapped container
   * `collRef`.
   * For example, a `collRef` of type `std::reference_wrapper` will result
   * in a C reference of the wrapped container, while a C pointer is left
   * unchanged and a `std::unique_ptr` object is turned into the equivalent
   * pointer to its elements.
   */
  template <typename CollRef>
  decltype(auto) collection_from_reference(CollRef& collRef);
  
  /// @}
  //--- END ContainerMetaprogramming -------------------------------------------
  
  
} // namespace util


//------------------------------------------------------------------------------
//--- template implementation
//------------------------------------------------------------------------------
namespace util {
  
  namespace details {
    
    //--------------------------------------------------------------------------
    //--- collection_value_XXXX
    //--------------------------------------------------------------------------
    template <typename Ptr>
    struct collection_value_type_impl_pointer;
    
    template <typename T>
    struct collection_value_type_impl_pointer<T*> {
      using type = T;
      using value_type = type;
    }; // struct collection_value_type_impl_pointer<T*>
    
    template <typename T>
    struct collection_value_type_impl_pointer<T[]> {
      using type = T;
      using value_type = type;
    }; // struct collection_value_type_impl_pointer<T[]>
    
    template <typename T, std::size_t N>
    struct collection_value_type_impl_pointer<T[N]> {
      using type = T;
      using value_type = type;
    }; // struct collection_value_type_impl_pointer<T[N]>
    
    
    template <typename Ptr, typename = void>
    struct collection_value_type_impl_unique_ptr
      : collection_value_type_impl_pointer<typename Ptr::pointer>
      {};
    
    
    template <typename Coll, typename = void>
    struct collection_value_type_impl {
      using type = typename Coll::value_type;
      using value_type = type;
    }; // struct collection_value_type_impl
    
    template <typename Coll>
    struct collection_value_type_impl
      <Coll, std::enable_if_t<std::is_pointer_v<std::decay_t<Coll>>>>
      : collection_value_type_impl_pointer<std::decay_t<Coll>>
    {};
    
    template <typename Coll>
    struct collection_value_type_impl
      <Coll, std::enable_if_t<util::is_unique_ptr_v<std::decay_t<Coll>>>>
      : collection_value_type_impl_unique_ptr<std::decay_t<Coll>>
    {};
    
    
    //--------------------------------------------------------------------------
    template <typename Coll, typename = void>
    struct collection_value_access_type_impl {
        private:
      static auto getBegin(Coll&& coll)
        { using std::begin; return begin(coll); }
      
        public:
      using type = decltype(*getBegin(std::declval<Coll>()));
      using value_type = collection_value_t<Coll>;
      
    }; // struct collection_value_access_type_impl
    
    template <typename T>
    struct collection_value_access_type_impl<T*, void> {
      using type = decltype(*(std::declval<T*>()));
      using value_type = T;
    }; // struct collection_value_access_type_impl<T*>
    
    
    template <typename Ptr>
    struct collection_value_access_type_impl
      <Ptr, std::enable_if_t<util::is_unique_ptr_v<std::decay_t<Ptr>>>>
     : collection_value_access_type_impl
       <std::remove_reference_t<typename Ptr::pointer>>
    {};
    
    
    //--------------------------------------------------------------------------
    template <typename Coll, typename = void>
    struct collection_value_constant_access_type_impl {
        private:
      static auto getCBegin(Coll&& coll)
        { using std::cbegin; return cbegin(coll); }
      
        public:
      using type = decltype(*getCBegin(std::declval<Coll>()));
      using value_type = collection_value_t<Coll>;
      
    }; // struct collection_value_constant_access_type_impl
    
    template <typename T>
    struct collection_value_constant_access_type_impl<T*, void> {
      using type = decltype(*(std::declval<std::add_const_t<T>*>()));
      using value_type = std::add_const_t<T>;
    }; // struct collection_value_constant_access_type_impl
    
    template <typename Ptr>
    struct collection_value_constant_access_type_impl
      <Ptr, std::enable_if_t<util::is_unique_ptr_v<std::decay_t<Ptr>>>>
     : collection_value_constant_access_type_impl
       <std::remove_reference_t<typename Ptr::pointer>>
    {};
    
    
    //--------------------------------------------------------------------------
    //--- util::make_collection_reference
    
    template <typename Coll, typename = void>
    struct make_collection_reference_impl {
      using type = std::reference_wrapper<Coll>;
      static auto make(Coll& coll) { return std::ref(coll); }
    }; // make_collection_reference_impl
    
    template <typename Coll>
    struct make_collection_reference_impl
      <Coll, std::enable_if_t<util::is_reference_wrapper_v<Coll>>>
    {
      using type = std::remove_cv_t<Coll>;
      static type make(Coll& refw) { return refw; }
    }; // make_collection_reference_impl<std::reference_wrapper>
    
    template <typename Coll>
    struct make_collection_reference_impl
      <Coll, std::enable_if_t<util::is_unique_ptr_v<Coll>>>
    {
      using type = typename Coll::pointer;
      static type make(Coll& uptr) { return uptr.get(); }
    }; // make_collection_reference_impl<std::unique_ptr>
    
    template <typename Ptr>
    struct make_collection_reference_impl
      <Ptr, std::enable_if_t<std::is_pointer_v<std::decay_t<Ptr>>>>
    {
      using type =
        std::add_pointer_t< // finally add the pointer
          std::remove_all_extents_t< // if it's a C array
            std::remove_pointer_t<std::decay_t<Ptr>>
          >
        >;
      static type make(Ptr& ptr) { return ptr; }
    }; // make_collection_reference_impl<std::unique_ptr>
    
    
    //--------------------------------------------------------------------------
    //--- util::collection_from_reference
    
    template <typename CollRef, typename = void>
    struct collection_from_reference_impl {
      using type
        = std::add_lvalue_reference_t<std::remove_reference_t<CollRef>>;
      static CollRef& get(CollRef& coll) { return coll; }
    }; // collection_from_reference_impl
    
    template <typename CollRef>
    struct collection_from_reference_impl
      <CollRef, std::enable_if_t<util::is_reference_wrapper_v<CollRef>>>
    {
      using type = std::add_lvalue_reference_t<typename CollRef::type>;
      static type get(CollRef& refw) { return refw.get(); }
    }; // collection_from_reference_impl<std::reference_wrapper>
    
    template <typename CollRef>
    struct collection_from_reference_impl
      <CollRef, std::enable_if_t<util::is_unique_ptr_v<CollRef>>>
    {
      using type = typename CollRef::pointer;
      static type get(CollRef& uptr) { return uptr.get(); }
    }; // collection_from_reference_impl<std::unique_ptr>
    
    template <typename T>
    struct collection_from_reference_impl<T*> {
      using type = T*;
      static type get(T* ptr) { return ptr; }
    }; // collection_from_reference_impl<T*>
    
    template <typename T>
    struct collection_from_reference_impl<T[]> {
      using type = T*;
      static type get(T ptr[]) { return ptr; }
    }; // collection_from_reference_impl<T[]>
    
    template <typename T, std::size_t N>
    struct collection_from_reference_impl<T[N]> {
      using type = T*;
      static type get(T ptr[N]) { return ptr; }
    }; // collection_from_reference_impl<T[N]>
    
    
    //--------------------------------------------------------------------------
    
  } // namespace details
  
  
  //----------------------------------------------------------------------------
  template <typename Coll>
  struct collection_value_type {
    // remove all referenceness, constantness etc. from `Coll`;
    // also remove all referenceness from the result
    using type = util::strip_referenceness_t<
      typename details::collection_value_type_impl
        <util::strip_referenceness_t<Coll>>::type
      >
      ;
  };
  
  
  //----------------------------------------------------------------------------
  template <typename Coll>
  struct collection_value_access_type
    : public details::collection_value_access_type_impl
      <util::strip_referenceness_t<Coll>>
  {};
  
  
  //----------------------------------------------------------------------------
  template <typename Coll>
  struct collection_value_constant_access_type
    : public details::collection_value_constant_access_type_impl
      <util::strip_referenceness_t<Coll>>
  {};
  
  
  //----------------------------------------------------------------------------
  template <typename Coll>
  struct collection_reference_type
    : details::make_collection_reference_impl<std::remove_reference_t<Coll>>
  {};

  //----------------------------------------------------------------------------
  template <typename Coll>
  auto make_collection_reference(Coll&& coll) { 
   return details::make_collection_reference_impl<std::remove_reference_t<Coll>>
     ::make(coll)
     ;
  }

  //----------------------------------------------------------------------------
  template <typename CollRef>
  struct collection_from_reference_type
    : details::collection_from_reference_impl<std::remove_reference_t<CollRef>>
    {};

  //----------------------------------------------------------------------------
  template <typename CollRef>
  decltype(auto) collection_from_reference(CollRef& collRef)
   { return details::collection_from_reference_impl<CollRef>::get(collRef); }

  //----------------------------------------------------------------------------
  
} // namespace util


//------------------------------------------------------------------------------

#endif // LARCOREALG_COREUTILS_CONTAINERMETA_H

