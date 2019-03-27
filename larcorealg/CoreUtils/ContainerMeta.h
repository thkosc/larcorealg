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
    template <typename Coll>
    struct collection_value_type_impl {
      using type = typename Coll::value_type;
      using value_type = type;
    }; // struct collection_value_type_impl
    
    template <typename T>
    struct collection_value_type_impl<T*> {
      using type = T;
      using value_type = type;
    }; // struct collection_value_type_impl<T*>
    
    template <typename T, std::size_t N>
    struct collection_value_type_impl<T[N]>: collection_value_type_impl<T*> {};
    
    
    //--------------------------------------------------------------------------
    template <typename Coll>
    struct collection_value_access_type_impl {
        private:
      static auto getBegin(Coll&& coll)
        { using std::begin; return begin(coll); }
      
        public:
      using type = decltype(*getBegin(std::declval<Coll>()));
      using value_type = collection_value_t<Coll>;
      
    }; // struct collection_value_access_type_impl
    
    
    //--------------------------------------------------------------------------
    template <typename Coll>
    struct collection_value_constant_access_type_impl {
        private:
      static auto getCBegin(Coll&& coll)
        { using std::cbegin; return cbegin(coll); }
      
        public:
      using type = decltype(*getCBegin(std::declval<Coll>()));
      using value_type = collection_value_t<Coll>;
      
    }; // struct collection_value_constant_access_type_impl
    
    
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
  
} // namespace util


//------------------------------------------------------------------------------

#endif // LARCOREALG_COREUTILS_CONTAINERMETA_H

