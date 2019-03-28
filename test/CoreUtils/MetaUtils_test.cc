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

BOOST_AUTO_TEST_CASE(NothingTestCase) {
} // BOOST_AUTO_TEST_CASE(NothingTestCase)
