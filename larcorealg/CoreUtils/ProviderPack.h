/**
 * @file    larcorealg/CoreUtils/ProviderPack.h
 * @brief   Data structure containing constant pointers to classes
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    November 20th, 2015
 * @version 1.0
 * 
 * This is a header-only library.
 * It depends only on standard C++ and does not require additional linkage.
 */


#ifndef LARCOREALG_COREUTILS_PROVIDERPACK_H
#define LARCOREALG_COREUTILS_PROVIDERPACK_H


// C/C++ standard library
#include <tuple>
#include <type_traits>
#include <limits>


namespace lar {
  
  namespace details {
    
    template <typename... Types>
    struct has_duplicate_types;
    
    
    /**
     * @brief Index of the class among Bases which is base of Derived.
     * @tparam Derived the class to be found
     * @tparam Bases a list of classes candidate to be the base of Derived
     * @return index of the class among Bases which is base of Derived
     * @throw static_assert if multiple classes are base of `Derived`
     * @see hasBaseOf(), findBaseOf()
     * 
     * If no class among `Bases` is actually a base class of `Derived`, an
     * invalid index is returned, greater than any valid index (that is,
     * no smaller than `sizeof...(Bases)`).
     */
    template <typename Derived, typename... Bases>
    constexpr std::size_t indexOfBaseOf();
    
    template <typename Derived, typename... Bases>
    constexpr std::size_t indexOfDerivedFrom();
    
    /**
     * @brief Index of the class among Bases which is base of Derived.
     * @tparam Derived the class to be found
     * @tparam Bases a list of classes candidate to be the base of Derived
     * @return index of the class among Bases which is base of Derived
     * @throw static_assert if none, or multiple classes, are base of `Derived`
     * @see hasBaseOf(), indexOfBaseOf()
     */
    template <typename Derived, typename... Bases>
    constexpr std::size_t findBaseOf();
    
    template <typename Derived, typename... Bases>
    constexpr std::size_t findDerivedFrom();
    
    /**
     * @brief Returns whether there is exactly one base class of `Derived` among
     *        `Bases`.
     * @tparam Derived the class to be found
     * @tparam Bases a list of classes candidate to be the base of Derived
     * @return whether there is exactly one base class of `Derived`
     * @throw static_assert if multiple classes are base of `Derived`
     * @see indexOfBaseOf(), findBaseOf()
     */
    template <typename Derived, typename... Bases>
    constexpr std::size_t hasBaseOf()
      { return indexOfBaseOf<Derived, Bases...>() < sizeof...(Bases); }
    
    template <typename Derived, typename... Bases>
    constexpr std::size_t hasDerivedFrom()
      { return indexOfDerivedFrom<Derived, Bases...>() < sizeof...(Bases); }
    
    
    /// Implementation detail for the extraction constructor
    template
      <typename DestPack, typename SourcePack, typename... ExtractProviders>
    struct SetFrom;
 
  } // namespace details
  
  
  /** **************************************************************************
   * @brief Container for a list of pointers to providers
   * @tparam Providers types of the providers in the parameter pack
   * 
   * The pointers are stored as constant.
   * Note that this container can host any type of objects, and it has 
   * "provider" in the name because the reason it was written was to provide
   * a fast way to specify a set of LArSoft service providers.
   * The only limitation is that there should be only one object per type.
   * Pointed objects are not owned by this class.
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
   * A a;
   * B1 b; // derived from B
   * C c;
   * D d;
   * ProviderPack<A, B, C> pack(&a, &b, &c);
   * 
   * // obtain a constant pointer to b from pack:
   * B const* b_ptr = pack.get<B>();
   * 
   * if constexpr (pack.has<D>()) std::cerr << "Unexpected!" << std::endl;
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * (note that in the latter check `constexpr` is supported only since C++17).
   */
  template <typename... Providers>
  class ProviderPack {
    static_assert(!details::has_duplicate_types<Providers...>::value,
      "Providers in ProviderPack are repeated");
    
    using this_type = ProviderPack<Providers...>; ///< alias of this class
    
    /// type used for storage of the pointers
    using tuple_type = std::tuple<Providers const*...>;
    
      public:
    
    /// Default constructor: a null provider pointer for each type
    ProviderPack() = default;
    
    /// Constructor: stores a provider pointer for each type
    ProviderPack(Providers const* ...provider_ptrs): providers(provider_ptrs...)
      {}
    
    /**
     * @brief Constructor: extracts the providers from another parameter pack
     * @tparam OtherProviders list of the providers of the source provider pack
     * @param from where to copy the information from
     * 
     * This constructor requires all the providers we need to be present
     * in the source provider pack.
     */
    template<typename... OtherProviders>
    ProviderPack(ProviderPack<OtherProviders...> const& from)
      {
        details::SetFrom
          <this_type, ProviderPack<OtherProviders...>, Providers...>
          (*this, from);
      }

    /**
     * @brief Constructor: picks the providers from the specified ones
     * @tparam OtherProviders list of the type of providers offered
     * @param providers all the providers needed (or more)
     * 
     * This constructor will pick, among the offered providers, the ones that
     * are needed.
     */
    template<typename... OtherProviders>
    ProviderPack(OtherProviders const*... providers)
      {
        details::SetFrom
          <this_type, ProviderPack<OtherProviders...>, Providers...>
          (*this, ProviderPack<OtherProviders...>(providers...));
      }

    /**
     * @brief Constructor: picks the providers from a pack plus specified ones
     * @tparam FromPack parameter pack to start from
     * @tparam OtherProviders list of the type of providers offered
     * @param fromPack providers to be picked
     * @param providers all the remaining providers needed (or more)
     * @see expandProviderPack()
     * 
     * This constructor will pick all the providers from the specified pack,
     * and the ones from the other providers.
     * This constructor can be used to "expand" from another provider:
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     *     A a;
     *     B b;
     *     C c;
     *     D d;
     *     ProviderPack<A, D> pack(&a, &d);
     *     ProviderPack<A, B, C, D> largerPack(pack, &c, &b);
     * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     */
    template<typename... PackProviders, typename... OtherProviders>
    ProviderPack(
      ProviderPack<PackProviders...> const& fromPack,
      OtherProviders const*... providers
      );
    
    
    /// Returns the provider with the specified type
    template <typename Provider>
    Provider const* get() const
      {
        constexpr auto providerIndex
          = details::findDerivedFrom<Provider, Providers...>();
        return std::get<providerIndex>(providers);
      } // get<>()
    
    
    /// Sets the provider with the specified type
    template <typename Provider>
    void set(Provider const* provider_ptr)
      {
        constexpr auto providerIndex
          = details::findDerivedFrom<Provider, Providers...>();
        std::get<providerIndex>(providers) = provider_ptr;
      } // set<>()
    
    /// Returns whether there is a provider with the specified type
    template <typename Provider>
    static constexpr bool has()
      { return details::hasDerivedFrom<Provider, Providers...>(); }
    
    
    /// Returns whether other provider pack has all the same providers as this
    template <typename... OtherProviders>
    bool operator== (ProviderPack<OtherProviders...> const& other) const;
    
    /// Returns whether other provider pack and this have different providers
    template <typename... OtherProviders>
    bool operator!= (ProviderPack<OtherProviders...> const& other) const;
    
    
    /**
     * @brief Returns whether all our providers are in the OfferedProviders list
     * @tparam OfferedProviders list of offered providers
     * 
     * This static function returns true if all the providers in this provider
     * pack are included among the OfferedProviders list. That list can contain
     * additional provider types, which will not affect the result.
     * 
     * Usage example:
     *     
     *     using providers_t
     *       = lar::ProviderPack<geo::GeometryCore, detinfo::LArProperties>;
     *     static_assert(
     *       providers_t::containsProviders
     *         <detinfo::LArProperties, detinfo::DetectorProperties>(),
     *       "Not all the required providers are present."
     *       );
     *     
     * In this example, the assertion will fail because of the absence of
     * `detinfo::DetectorProperties` from providers_t.
     */
    template <typename... OtherProviders>
    static constexpr bool containsProviders();
    
      private:
    
    tuple_type providers; ///< container of the pointers, type-safe
    
  }; // class ProviderPack
  
  
  /**
   * @brief Function to create a ProviderPack from the function arguments
   * @tparam Providers types of the providers in the parameter pack
   * @param providers constant pointers to the providers
   * @return a ProviderPack object containing all the specified providers
   *
   * This is an convenience function to reduce the typing needed to instantiate
   * a ProviderPack. Example:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *     A a;
   *     B b;
   *     auto pack = makeProviderPack(&a, &b);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * creates a `ProviderPack<A, B>`.
   */
  template <typename... Providers>
  ProviderPack<Providers...> makeProviderPack(Providers const* ...providers)
    { return ProviderPack<Providers...>(providers...); }
  
  
  /**
   * @brief Function to create a ProviderPack by adding to another
   * @tparam PackProviders types of the providers in the original parameter pack
   * @tparam MoreProviders types of the providers to be added
   * @param pack parameter pack with the first providers
   * @param providers constant pointers to the other providers to be added
   * @return a ProviderPack object containing all the specified providers
   *
   * This is an convenience function to reduce the typing needed to instantiate
   * a ProviderPack. Use it like:
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   *     A a;
   *     B b;
   *     C c;
   *     D d;
   *     auto pack = makeProviderPack(&a, &d);
   *     auto largerPack = expandProviderPack(pack, &c, &b);
   * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   * creates a `ProviderPack<A, D, C, B>` including all the four objects.
   */
  template <typename... PackProviders, typename... MoreProviders>
  ProviderPack<PackProviders..., MoreProviders...> expandProviderPack(
    ProviderPack<PackProviders...> const& pack,
    MoreProviders const* ...providers
    )
    { return { pack, providers... }; }
  
  
} // namespace lar


//------------------------------------------------------------------------------
//--- Implementation details
//---
namespace lar {
  namespace details {
    
    template <bool Value>
    using bool_constant = std::integral_constant<bool, Value>; // also in C++17
    
    template <std::size_t Value>
    using index_constant = std::integral_constant<std::size_t, Value>;
    
    //--------------------------------------------------------------------------
    //--- has_type, has_duplicate_types, are_same_types
    //---
    template <typename Target, typename... Types>
    struct has_type;
    
    template <typename Target, typename First, typename... Others>
    struct has_type<Target, First, Others...>: has_type<Target, Others...> {};
    
    template <typename Target, typename... Others>
    struct has_type<Target, Target, Others...>: std::true_type {};
    
    template <typename Target>
    struct has_type<Target>: std::false_type {};
    
    
    //--------------------------------------------------------------------------
    template <typename Key, typename... Types>
    struct has_duplicate_types<Key, Types...>:
      public bool_constant
        <has_type<Key, Types...>() || has_duplicate_types<Types...>()>
      {};
    
    template <>
    struct has_duplicate_types<>: public std::false_type {};
    
    
    //--------------------------------------------------------------------------
    template <typename... Types>
    struct are_types_contained;
    
    template <typename... Types>
    struct are_same_types {
      
      template <typename... AsTypes>
      static constexpr bool as()
        {
          return (sizeof...(Types) == sizeof...(AsTypes))
            && are_types_contained<Types...>::template in<AsTypes...>();
        }
      
    }; // are_same_types
    
    
    template <typename First, typename... OtherTypes>
    struct are_types_contained<First, OtherTypes...> {
      template <typename... AsTypes>
      static constexpr bool in()
        {
          return are_types_contained<OtherTypes...>::template in<AsTypes...>()
            && has_type<First, AsTypes...>();
        }
    };
    
    template <typename First>
    struct are_types_contained<First> {
      template <typename... AsTypes>
      static constexpr bool in()
        { return has_type<First, AsTypes...>(); }
    };

    
    template <typename T>
    struct is_provider_pack: public std::false_type {};
    
    template <typename... Providers>
    struct is_provider_pack<ProviderPack<Providers...>>: public std::true_type
      {};
    
    
    template <typename APack, typename BPack>
    struct have_same_provider_types: public std::false_type {};
    
    template <typename... AProviders, typename... BProviders>
    struct have_same_provider_types
      <ProviderPack<AProviders...>, ProviderPack<BProviders...>>
      : public std::integral_constant
        <bool, are_same_types<AProviders...>::template as<BProviders...>()>
      {};
    
    
    // --- BEGIN impementation of findDerivedFrom() ----------------------------
    // 
    // This is the implementation of findDerivedFrom() and findBaseOf().
    // The functions are expected to assert that there is exactly one answer to
    // a matching condition: for findDerivedFrom() that answer is which is the
    // derived class of Base among Derived, for findBaseOf() is the opposite.
    // 
    // The implementation finds and returns the index of first class matching
    // the condition, and it asserts that the answer is valid.
    // It also finds the index of the next class matching the condition, just to
    // assert that it is not valid (that is, there is no other matching class).
    // 
    // The class returning the index of the first matching class is implemented
    // recursively, on a class taking the target class, the first of the
    // candidate classes, and then all the others in a parameter pack.
    // When the first candidate satisfies the condition, then the recursion is
    // over. In the simplest case, where we are just looking for a specific
    // class, the condition (the Target is the same as the First Candidate) can
    // be implied and the compiler can implement it directly with parameter
    // matching. In our cases the condition is non-trivial and we use a
    // two-level recursion: the outer level, the "dispatcher", invokes a
    // different "answer" class whether the condition is true or false. In the
    // former case a result is presented; in the other case recursion ensues,
    // back to the dispatcher. The dispatcher can also terminate the recursion
    // when no answer is found. If recursion is terminated without a match, an
    // invalid index is returned.
    // 
    // The class returning the next matching class is simply skipping a number
    // of candidates, and then behaving like the one looking for the first
    // matching candidate.
    //
    
    //
    // class to find the first matching class
    //
    template <
      template <typename A, typename B> class Match,
      typename Target, bool IsMatch, typename... Candidates
      >
    struct findFirstMatching_answer;
    
    template <
      template <typename A, typename B> class Match,
      typename Target, typename... Candidates
      >
    struct findFirstMatching_answer<Match, Target, true, Candidates...>
      : public index_constant<0U>
    {};
    
    template <
      template <typename A, typename B> class Match,
      typename Target, typename... Candidates
      >
    struct findFirstMatching_dispatcher;
    
    // end-of-recursion
    template <
      template <typename A, typename B> class Match,
      typename Target
      >
    struct findFirstMatching_dispatcher<Match, Target>
      : findFirstMatching_answer<Match, Target, true>
    {};
    
    template <
      template <typename A, typename B> class Match,
      typename Target, typename FirstCandidate, typename... OtherCandidates
      >
    struct findFirstMatching_dispatcher
      <Match, Target, FirstCandidate, OtherCandidates...>
      : public findFirstMatching_answer<
        Match,
        Target,
        Match<FirstCandidate, Target>::value,
        FirstCandidate,
        OtherCandidates...
      >
    {};
    
    template <
      template <typename A, typename B> class Match,
      typename Target, typename FirstCandidate, typename... OtherCandidates
      >
    struct findFirstMatching_answer
      <Match, Target, false, FirstCandidate, OtherCandidates...>
      : public index_constant
        <(1U + findFirstMatching_dispatcher<Match, Target, OtherCandidates...>::value)>
    {};
    
    template <
      template <typename A, typename B> class Match,
      typename Target, typename... Candidates
      >
    struct findFirstMatching_impl
      : findFirstMatching_dispatcher<Match, Target, Candidates...>
    {
        private:
      static constexpr auto _index
        = findFirstMatching_dispatcher<Match, Target, Candidates...>();
    }; // struct findFirstMatching_impl
    
    
    //
    // class to apply findFirstMatching_impl after skipping some candidates
    //
    template <
      unsigned int NSkip,
      template <typename A, typename B> class Match,
      typename Target, typename... Candidates
      >
    struct findNextMatching_impl;
    
    // recursion: peel one
    template <
      unsigned int NSkip,
      template <typename A, typename B> class Match,
      typename Target, typename FirstCandidate, typename... OtherCandidates
      >
    struct findNextMatching_impl
      <NSkip, Match, Target, FirstCandidate, OtherCandidates...>
      : index_constant<(
        1U
        + findNextMatching_impl
          <(NSkip - 1U), Match, Target, OtherCandidates...>::value
      )>
    {
      static_assert(NSkip > 0U, "Implementation error: no arguments to skip!");
    };
    
    // end-of-recursion: skipped enough
    template <
      template <typename A, typename B> class Match,
      typename Target, typename FirstCandidate, typename... OtherCandidates
      >
    struct findNextMatching_impl
      <0U, Match, Target, FirstCandidate, OtherCandidates...>
      : findFirstMatching_impl
        <Match, Target, FirstCandidate, OtherCandidates...>
    {};
    
    // end-of-recursion: all arguments skipped
    template <
      unsigned int NSkip,
      template <typename A, typename B> class Match,
      typename Target
      >
    struct findNextMatching_impl<NSkip, Match, Target>
      : findFirstMatching_impl<Match, Target>
      {};
    
    //
    // class finding a match and asserting its existence and unicity
    //
    template <
      template <typename A, typename B> class Match,
      typename Target, typename... Candidates
      >
    struct findTheMatching_impl
      : findFirstMatching_impl<Match, Target, Candidates...>
    {
        private:
      static constexpr auto _index
        = findFirstMatching_dispatcher<Match, Target, Candidates...>();
      
      static_assert(
        findNextMatching_impl<_index + 1U, Match, Target, Candidates...>()
          >= sizeof...(Candidates),
        "Multiple candidate classes match the Target one"
        );
    }; // struct findTheMatching_impl
    
    //
    // implementations with concrete matching conditions
    //
    template <typename Derived, typename... Bases>
    constexpr std::size_t indexOfBaseOf()
      { return findTheMatching_impl<std::is_base_of, Derived, Bases...>(); }
    
    template <typename Derived, typename... Bases>
    constexpr std::size_t findBaseOf()
      { 
        constexpr std::size_t index = indexOfBaseOf<Derived, Bases...>(); 
        static_assert(
          index < sizeof...(Bases),
          "Target is not derived from any of the available classes"
          );
        return index;
      } // findBaseOf()
    
    // this matching condition is the mirror of std::is_base_of
    template <typename Derived, typename Base>
    struct is_derived_of: std::is_base_of<Base, Derived> {};
    
    template <typename Base, typename... Derived>
    constexpr std::size_t indexOfDerivedFrom()
      { return findTheMatching_impl<is_derived_of, Base, Derived...>(); }
    
    template <typename Base, typename... Derived>
    constexpr std::size_t findDerivedFrom()
      { 
        constexpr std::size_t index = indexOfDerivedFrom<Base, Derived...>(); 
        static_assert(
          index < sizeof...(Derived),
          "Target is not base of any of the available classes"
          );
        return index;
      } // findDerivedFrom()
    
    // --- END impementation of findDerivedFrom() ------------------------------
    
    
    //--------------------------------------------------------------------------
    //--- SetFrom
    //---
    template <
      typename DestPack, typename SourcePack,
      typename FirstProvider, typename... OtherProviders
      >
    struct SetFrom<DestPack, SourcePack, FirstProvider, OtherProviders...> {
      SetFrom(DestPack& pack, SourcePack const& from)
        {
          pack.set(from.template get<FirstProvider>());
          SetFrom<DestPack, SourcePack, OtherProviders...>(pack, from);
        }
    }; // SetFrom<First, Others...>
    
    template <typename DestPack, typename SourcePack>
    struct SetFrom<DestPack, SourcePack> {
      SetFrom(DestPack&, SourcePack const&) {}
    };

    //--------------------------------------------------------------------------
    //--- Compare
    //---
    template <typename Provider, typename APack, typename BPack>
    bool haveSameProvider(APack const& a, BPack const& b) {
      static_assert(is_provider_pack<APack>() && is_provider_pack<BPack>(),
        "This class needs two ProviderPack template types.");
      return a.template get<Provider>() == b.template get<Provider>();
    } // haveSameProvider()
    
    
    template <typename APack, typename BPack>
    struct ProviderPackComparerBase {
      
      static_assert(have_same_provider_types<APack, BPack>(),
        "The specified provider packs have different types.");
      
    }; // ProviderPackComparerBase
    
    
    template <typename APack, typename BPack, typename... Providers>
    struct ProviderPackComparer;
    
    template
      <typename APack, typename BPack, typename First, typename... Others>
    struct ProviderPackComparer<APack, BPack, First, Others...>
      : ProviderPackComparerBase<APack, BPack>
    {
      static bool compare (APack const& a, BPack const& b)
        {
          return haveSameProvider<First>(a, b)
            && ProviderPackComparer<APack, BPack, Others...>::compare(a, b);
        }
    }; // ProviderPackComparer<APack, BPack, First, Others...>
    
    template
      <typename APack, typename BPack, typename First>
    struct ProviderPackComparer<APack, BPack, First>
      : ProviderPackComparerBase<APack, BPack>
    {
      static bool compare (APack const& a, BPack const& b)
        { return haveSameProvider<First>(a, b); }
    }; // ProviderPackComparer<APack, BPack, First>
    
    
    //--------------------------------------------------------------------------

  } // namespace details
  
  
  //----------------------------------------------------------------------------
  //--- ProviderPack
  //---
  
  template <typename... Providers>
  template<typename... PackProviders, typename... OtherProviders>
  ProviderPack<Providers...>::ProviderPack(
    ProviderPack<PackProviders...> const& fromPack,
    OtherProviders const*... providers
    )
  {
    
    // verify that the list of providers in argument is the exact one we need
    static_assert(
      details::are_same_types<Providers...>
        ::template as<PackProviders..., OtherProviders...>(),
      "The providers types in the arguments do not match the ones needed."
      );
    
    // copy all the providers from the provider pack
    details::SetFrom
      <this_type, ProviderPack<PackProviders...>, PackProviders...>
      (*this, fromPack);
    
    // put the other providers in a temporary parameter pack, and copy it
    // (this is convenience, a direct implementation would be probably better)
    details::SetFrom
      <this_type, ProviderPack<OtherProviders...>, OtherProviders...>
      (*this, makeProviderPack(providers...));
    
  } // ProviderPack<Providers...>::ProviderPack(ProviderPack, OtherProviders...)
  
  
  //----------------------------------------------------------------------------
  template <typename... Providers>
  template <typename... OtherProviders>
  bool ProviderPack<Providers...>::operator==
    (ProviderPack<OtherProviders...> const& other) const
  {
    return details::ProviderPackComparer<
      ProviderPack<Providers...>, ProviderPack<OtherProviders...>, Providers...
      >::compare(*this, other);
  }
  
  
  template <typename... Providers>
  template <typename... OtherProviders>
  bool ProviderPack<Providers...>::operator!=
    (ProviderPack<OtherProviders...> const& other) const
    { return !(*this == other); }
  
 
  //----------------------------------------------------------------------------
  template <typename... Providers>
  template <typename... OfferedProviders>
  constexpr bool ProviderPack<Providers...>::containsProviders() {
    return details::are_types_contained<Providers...>
      ::template in<OfferedProviders...>();
  } // ProviderPack<>::containsProviders()
  
  
  //----------------------------------------------------------------------------
  
} // namespace lar

#endif // LARCOREALG_COREUTILS_PROVIDERPACK_H
