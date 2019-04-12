/**
 * @file    MetaUtils_test.cc
 * @brief   Unit test for some of the utilities in MetaUtils.h
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    March 22, 2019
 * @see     `larcorealg/CoreUtils/MetaUtils.h`
 */

// Boost libraries
#define BOOST_TEST_MODULE ( MetaUtils_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// LArSoft libraries
#include "larcorealg/CoreUtils/MetaUtils.h"

// C/C++ standard libraries
#include <array>
#include <memory> // std::unique_ptr, std::shared_ptr
#include <type_traits>


//------------------------------------------------------------------------------
//--- static tests (would fail at compile time)
//------------------------------------------------------------------------------

static_assert( util::is_instance_of_v<std::unique_ptr, std::unique_ptr<int>>);
static_assert( util::is_instance_of_v<std::unique_ptr, std::unique_ptr<int> const&>);
static_assert(!util::is_instance_of_v<std::unique_ptr, std::shared_ptr<int>>);
static_assert(!util::is_instance_of_v<std::unique_ptr, int>);

//------------------------------------------------------------------------------
// util::is_STLarray_v
//------------------------------------------------------------------------------

static_assert(!util::is_STLarray_v<int>);
static_assert(!util::is_STLarray_v<int[3]>);
static_assert( util::is_STLarray_v<std::array<int, 3>>);


//------------------------------------------------------------------------------
// util::is_reference_wrapper_v
//------------------------------------------------------------------------------

static_assert(!util::is_reference_wrapper_v<int>);
static_assert( util::is_reference_wrapper_v<std::reference_wrapper<int>>);
static_assert( util::is_reference_wrapper_v<std::reference_wrapper<int> const&>);


//------------------------------------------------------------------------------
// util::is_unique_ptr_v
//------------------------------------------------------------------------------

static_assert(!util::is_unique_ptr_v<int>);
static_assert(!util::is_unique_ptr_v<int*>);
static_assert(!util::is_unique_ptr_v<int const*>);
static_assert( util::is_unique_ptr_v<std::unique_ptr<int>>);
static_assert( util::is_unique_ptr_v<std::unique_ptr<int> const>);
static_assert( util::is_unique_ptr_v<std::unique_ptr<int> const&>);
static_assert( util::is_unique_ptr_v<std::unique_ptr<int const>&>);


//------------------------------------------------------------------------------
// util::with_const_as_t
//------------------------------------------------------------------------------

// simple type
static_assert(std::is_same_v<util::with_const_as_t<double               , int      >, double               >);
static_assert(std::is_same_v<util::with_const_as_t<double          const, int      >, double               >);
static_assert(std::is_same_v<util::with_const_as_t<double               , int const>, double          const>);
static_assert(std::is_same_v<util::with_const_as_t<double          const, int const>, double          const>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile      , int      >, double volatile      >);
static_assert(std::is_same_v<util::with_const_as_t<double volatile const, int      >, double volatile      >);
static_assert(std::is_same_v<util::with_const_as_t<double volatile      , int const>, double volatile const>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile const, int const>, double volatile const>);

// l-value reference type
static_assert(std::is_same_v<util::with_const_as_t<double               &, int      >, double               &>);
static_assert(std::is_same_v<util::with_const_as_t<double          const&, int      >, double               &>);
static_assert(std::is_same_v<util::with_const_as_t<double               &, int const>, double          const&>);
static_assert(std::is_same_v<util::with_const_as_t<double          const&, int const>, double          const&>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile      &, int      >, double volatile      &>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile const&, int      >, double volatile      &>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile      &, int const>, double volatile const&>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile const&, int const>, double volatile const&>);

// r-value reference type
static_assert(std::is_same_v<util::with_const_as_t<double               &&, int      >, double               &&>);
static_assert(std::is_same_v<util::with_const_as_t<double          const&&, int      >, double               &&>);
static_assert(std::is_same_v<util::with_const_as_t<double               &&, int const>, double          const&&>);
static_assert(std::is_same_v<util::with_const_as_t<double          const&&, int const>, double          const&&>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile      &&, int      >, double volatile      &&>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile const&&, int      >, double volatile      &&>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile      &&, int const>, double volatile const&&>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile const&&, int const>, double volatile const&&>);


// simple type (key is reference)
static_assert(std::is_same_v<util::with_const_as_t<double               , int      &>, double               >);
static_assert(std::is_same_v<util::with_const_as_t<double          const, int      &>, double               >);
static_assert(std::is_same_v<util::with_const_as_t<double               , int const&>, double          const>);
static_assert(std::is_same_v<util::with_const_as_t<double          const, int const&>, double          const>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile      , int      &>, double volatile      >);
static_assert(std::is_same_v<util::with_const_as_t<double volatile const, int      &>, double volatile      >);
static_assert(std::is_same_v<util::with_const_as_t<double volatile      , int const&>, double volatile const>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile const, int const&>, double volatile const>);

// l-value reference type
static_assert(std::is_same_v<util::with_const_as_t<double               &, int      &>, double               &>);
static_assert(std::is_same_v<util::with_const_as_t<double          const&, int      &>, double               &>);
static_assert(std::is_same_v<util::with_const_as_t<double               &, int const&>, double          const&>);
static_assert(std::is_same_v<util::with_const_as_t<double          const&, int const&>, double          const&>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile      &, int      &>, double volatile      &>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile const&, int      &>, double volatile      &>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile      &, int const&>, double volatile const&>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile const&, int const&>, double volatile const&>);

// r-value reference type
static_assert(std::is_same_v<util::with_const_as_t<double               &&, int      &>, double               &&>);
static_assert(std::is_same_v<util::with_const_as_t<double          const&&, int      &>, double               &&>);
static_assert(std::is_same_v<util::with_const_as_t<double               &&, int const&>, double          const&&>);
static_assert(std::is_same_v<util::with_const_as_t<double          const&&, int const&>, double          const&&>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile      &&, int      &>, double volatile      &&>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile const&&, int      &>, double volatile      &&>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile      &&, int const&>, double volatile const&&>);
static_assert(std::is_same_v<util::with_const_as_t<double volatile const&&, int const&>, double volatile const&&>);

//------------------------------------------------------------------------------
//--- util::strip_referenceness_type
//------------------------------------------------------------------------------
static_assert(std::is_same_v<util::strip_referenceness_t<int                                                        >, int      >);
static_assert(std::is_same_v<util::strip_referenceness_t<int                                                 const  >, int const>);
static_assert(std::is_same_v<util::strip_referenceness_t<int                                                 const& >, int const>);
static_assert(std::is_same_v<util::strip_referenceness_t<int                                                      &&>, int      >);
static_assert(std::is_same_v<util::strip_referenceness_t<std::reference_wrapper<int                        >        >, int      >);
static_assert(std::is_same_v<util::strip_referenceness_t<std::reference_wrapper<int                        > const& >, int      >);
static_assert(std::is_same_v<util::strip_referenceness_t<std::reference_wrapper<std::reference_wrapper<int>>        >, int      >);


//------------------------------------------------------------------------------
//--- util::lvalue_reference_into_wrapper_type
//------------------------------------------------------------------------------
static_assert(std::is_same_v<util::lvalue_reference_into_wrapper_t<                       int               >,                        int             >);
static_assert(std::is_same_v<util::lvalue_reference_into_wrapper_t<                       int const         >,                        int const       >);
static_assert(std::is_same_v<util::lvalue_reference_into_wrapper_t<                       int             & >, std::reference_wrapper<int      >      >);
static_assert(std::is_same_v<util::lvalue_reference_into_wrapper_t<                       int const       & >, std::reference_wrapper<int const>      >);
static_assert(std::is_same_v<util::lvalue_reference_into_wrapper_t<                       int             &&>,                        int             >);
static_assert(std::is_same_v<util::lvalue_reference_into_wrapper_t<                       int const       &&>,                        int const       >);
static_assert(std::is_same_v<util::lvalue_reference_into_wrapper_t<std::reference_wrapper<int      >        >, std::reference_wrapper<int      >      >);
static_assert(std::is_same_v<util::lvalue_reference_into_wrapper_t<std::reference_wrapper<int      > const  >, std::reference_wrapper<int      > const>);
static_assert(std::is_same_v<util::lvalue_reference_into_wrapper_t<std::reference_wrapper<int const>        >, std::reference_wrapper<int const>      >);
static_assert(std::is_same_v<util::lvalue_reference_into_wrapper_t<std::reference_wrapper<int const> const  >, std::reference_wrapper<int const> const>);
static_assert(std::is_same_v<util::lvalue_reference_into_wrapper_t<std::reference_wrapper<int      >      & >, std::reference_wrapper<int      >      >);
static_assert(std::is_same_v<util::lvalue_reference_into_wrapper_t<std::reference_wrapper<int      > const& >, std::reference_wrapper<int      > const>);
static_assert(std::is_same_v<util::lvalue_reference_into_wrapper_t<std::reference_wrapper<int const>      & >, std::reference_wrapper<int const>      >);
static_assert(std::is_same_v<util::lvalue_reference_into_wrapper_t<std::reference_wrapper<int const> const& >, std::reference_wrapper<int const> const>);


//------------------------------------------------------------------------------
void referenced_address_test() {

  int v;
  int&       ref   = v;
  int const& cref  = v;
  auto       refw  = std::ref(v);
  auto       crefw = std::cref(v);

  BOOST_CHECK_EQUAL(util::referenced_address(ref  ), std::addressof(v));
  BOOST_CHECK_EQUAL(util::referenced_address(cref ), std::addressof(v));
  BOOST_CHECK_EQUAL(util::referenced_address(refw ), std::addressof(v));
  BOOST_CHECK_EQUAL(util::referenced_address(crefw), std::addressof(v));

} // lvalue_reference_into_wrapper_test()


//------------------------------------------------------------------------------
void lvalue_reference_into_wrapper_test() {

  int        obj   = 1;
  int      & ref   = obj;
  int const& cref  = obj;
  auto       refw  = std::ref(obj);
  auto       crefw = std::cref(obj);

  decltype(auto) ref_obj   = util::lvalue_reference_into_wrapper(obj  );
  decltype(auto) ref_ref   = util::lvalue_reference_into_wrapper(ref  );
  decltype(auto) ref_cref  = util::lvalue_reference_into_wrapper(cref );
  decltype(auto) ref_refw  = util::lvalue_reference_into_wrapper(refw );
  decltype(auto) ref_crefw = util::lvalue_reference_into_wrapper(crefw);

  // since we have already verified the correctness of
  // util::lvalue_reference_into_wrapper_t<>, we just check for consistency:
  static_assert(std::is_same_v<decltype(ref_obj  ), util::lvalue_reference_into_wrapper_t<int&           >>);
  static_assert(std::is_same_v<decltype(ref_ref  ), util::lvalue_reference_into_wrapper_t<decltype(ref  )>>);
  static_assert(std::is_same_v<decltype(ref_cref ), util::lvalue_reference_into_wrapper_t<decltype(cref )>>);
  static_assert(std::is_same_v<decltype(ref_refw ), util::lvalue_reference_into_wrapper_t<decltype(refw )>>);
  static_assert(std::is_same_v<decltype(ref_crefw), util::lvalue_reference_into_wrapper_t<decltype(crefw)>>);

  BOOST_CHECK_EQUAL(util::referenced_address(ref_obj  ), util::referenced_address(obj));
  BOOST_CHECK_EQUAL(util::referenced_address(ref_ref  ), util::referenced_address(obj));
  BOOST_CHECK_EQUAL(util::referenced_address(ref_cref ), util::referenced_address(obj));
  BOOST_CHECK_EQUAL(util::referenced_address(ref_refw ), util::referenced_address(obj));
  BOOST_CHECK_EQUAL(util::referenced_address(ref_crefw), util::referenced_address(obj));

} // lvalue_reference_into_wrapper_test()


//------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(ReferencesTestCase) {

  referenced_address_test();
  lvalue_reference_into_wrapper_test();

} // BOOST_AUTO_TEST_CASE(ReferencesTestCase)

//------------------------------------------------------------------------------
