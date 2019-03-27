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
#include <type_traits>
#include <array>


//------------------------------------------------------------------------------
//--- static tests (would fail at compile time)
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

BOOST_AUTO_TEST_CASE(NothingTestCase) {
} // BOOST_AUTO_TEST_CASE(NothingTestCase)
