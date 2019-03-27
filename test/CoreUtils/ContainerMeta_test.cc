/**
 * @file    ContainerMeta_test.cc
 * @brief   Unit test for some of the utilities in `ContainerMeta.h`
 * @author  Gianluca Petrillo (petrillo@fnal.gov)
 * @date    March 25, 2019
 * @see     `larcorealg/CoreUtils/ContainerMeta.h`
 */

// Boost libraries
#define BOOST_TEST_MODULE ( ContainerMeta_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// LArSoft libraries
#include "larcorealg/CoreUtils/ContainerMeta.h"

// C/C++ standard libraries
#include <list>
#include <functional> // std::reference_wrapper
#include <type_traits>


//------------------------------------------------------------------------------
//--- static tests (would fail at compile time)
//------------------------------------------------------------------------------
struct O {
  int value = 0;
};

static_assert(std::is_same_v<util::collection_value_t                <std::list<O                              >      >, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::list<O      *                       >      >, O      *       >);
static_assert(std::is_same_v<util::collection_value_t                <std::list<O const*                       >      >, O const*       >);

static_assert(std::is_same_v<util::collection_value_access_t         <std::list<O                              >      >, O             &>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::list<O      *                       >      >, O      *      &>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::list<O const*                       >      >, O const*      &>);

static_assert(std::is_same_v<util::collection_value_constant_access_t<std::list<O                              >      >, O        const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::list<O      *                       >      >, O      * const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::list<O const*                       >      >, O const* const&>);

static_assert(std::is_same_v<util::collection_value_t                <std::list<O                              > const>, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::list<O      *                       > const>, O      *       >);
static_assert(std::is_same_v<util::collection_value_t                <std::list<O const*                       > const>, O const*       >);

static_assert(std::is_same_v<util::collection_value_access_t         <std::list<O                              > const>, O        const&>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::list<O      *                       > const>, O      * const&>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::list<O const*                       > const>, O const* const&>);

static_assert(std::is_same_v<util::collection_value_constant_access_t<std::list<O                              > const>, O        const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::list<O      *                       > const>, O      * const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::list<O const*                       > const>, O const* const&>);

//
// std::reference_wrapper content
//
static_assert(std::is_same_v<util::collection_value_t                <std::list<std::reference_wrapper<O      >>      >, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::list<std::reference_wrapper<O      >> const>, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::list<std::reference_wrapper<O const>>      >, O const        >);
static_assert(std::is_same_v<util::collection_value_t                <std::list<std::reference_wrapper<O const>> const>, O const        >);

//
// std::reference_wrapper container
//
static_assert(std::is_same_v<util::collection_value_t                <                       std::list<O      >&      >, O              >);
static_assert(std::is_same_v<util::collection_value_t                <                       std::list<O      > const&>, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::reference_wrapper<std::list<O      >>      >, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::reference_wrapper<std::list<O      >>     &>, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::reference_wrapper<std::list<O      >>const&>, O              >);
static_assert(std::is_same_v<util::collection_value_t                <std::reference_wrapper<std::list<O> const>const&>, O              >);

static_assert(std::is_same_v<util::collection_value_access_t         <                       std::list<O      >&      >, O             &>);
static_assert(std::is_same_v<util::collection_value_access_t         <                       std::list<O      > const&>, O        const&>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::reference_wrapper<std::list<O      >>      >, O             &>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::reference_wrapper<std::list<O      >>     &>, O             &>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::reference_wrapper<std::list<O      >>const&>, O             &>);
static_assert(std::is_same_v<util::collection_value_access_t         <std::reference_wrapper<std::list<O> const>const&>, O        const&>);

static_assert(std::is_same_v<util::collection_value_constant_access_t<                       std::list<O      >&      >, O        const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<                       std::list<O      > const&>, O        const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::reference_wrapper<std::list<O      >>      >, O        const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::reference_wrapper<std::list<O      >>     &>, O        const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::reference_wrapper<std::list<O      >>const&>, O        const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::reference_wrapper<std::list<O> const>const&>, O        const&>);


//------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(NothingTestCase) {
} // BOOST_AUTO_TEST_CASE(NothingTestCase)
