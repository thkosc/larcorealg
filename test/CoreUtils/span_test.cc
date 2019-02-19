/**
 * @file   span_test.cc
 * @brief  Unit test for `util::span` class.
 * @author Gianluca Petrillo (petrillo@slac.stanford.edu)
 * @date   January 29, 2019
 * @see    `larcorealg/CoreUtils/span.h`
 */

// Boost libraries
#define BOOST_TEST_MODULE ( RealComparisons_test )
#include <cetlib/quiet_unit_test.hpp> // BOOST_AUTO_TEST_CASE()
#include <boost/test/test_tools.hpp> // BOOST_CHECK(), BOOST_CHECK_EQUAL()

// LArSoft libraries
#include "larcorealg/CoreUtils/span.h"

// C/C++ standard libraries
#include <iostream> // std::cout
#include <vector>
#include <numeric> // std::iota()
#include <type_traits> // std::is_same<>


//------------------------------------------------------------------------------
template <typename Iter, typename Cont>
void test_span(Cont& c) {
  using TestVector_t = Cont;
  
  TestVector_t v;
  
  auto r = util::make_span(v);
  
  // check the type of the object
  static_assert
    (std::is_same<decltype(r), util::span<Iter>>());
  static_assert(std::is_same<decltype(r), util::span<Iter, Iter>>());
  static_assert(
    std::is_same<typename decltype(r)::value_type, typename Iter::value_type>()
    );
  static_assert
    (std::is_same<typename decltype(r)::reference, typename Iter::reference>());
  
  BOOST_CHECK(r.begin() == v.begin());
  BOOST_CHECK(r.end() == v.end());
  
  BOOST_CHECK_EQUAL(r.empty(), v.empty());
  BOOST_CHECK_EQUAL(r.size(), v.size());
  
  auto iV = v.cbegin();
  unsigned int n = 0;
  for (auto i: r) {
    BOOST_CHECK_EQUAL(i, *iV);
    ++n;
    ++iV;
  } // for
  BOOST_CHECK_EQUAL(n, v.size());
  
} // test_span()


//------------------------------------------------------------------------------
template <typename Iter, typename Cont>
void test_const_span(Cont& c) {
  using TestVector_t = Cont;
  
  TestVector_t v;
  
  auto r = util::make_const_span(v);
  
  // check the type of the object
  static_assert
    (std::is_same<decltype(r), util::span<Iter>>());
  static_assert(std::is_same<decltype(r), util::span<Iter, Iter>>());
  static_assert(
    std::is_same<typename decltype(r)::value_type, typename Iter::value_type>()
    );
  static_assert
    (std::is_same<typename decltype(r)::reference, typename Iter::reference>());
  
  BOOST_CHECK(r.begin() == v.cbegin());
  BOOST_CHECK(r.end() == v.cend());
  
  BOOST_CHECK_EQUAL(r.empty(), v.empty());
  BOOST_CHECK_EQUAL(r.size(), v.size());
  
  auto iV = v.cbegin();
  unsigned int n = 0;
  for (auto i: r) {
    BOOST_CHECK_EQUAL(i, *iV);
    ++n;
    ++iV;
  } // for
  BOOST_CHECK_EQUAL(n, v.size());
  
} // test_const_span()


//------------------------------------------------------------------------------

struct RangeDocumentation1TestClass {
  using Data_t = std::vector<int>;

  int sum(util::span<Data_t::iterator> r) {
    decltype(r)::value_type s = 0;
    for (auto v: r) s += v;
    return s;
  } // sum

  void analyse() {
    
    Data_t v(10U);
    std::iota(v.begin(), v.end(), 1); // { 1, 2, 3, 4 ,5 ,6, 7, 8, 9, 10 }
    
    auto span5 = util::span(v.begin(), v.begin() + 5);
    std::cout << "Sum of 5 numbers: " << sum(span5) << std::endl;
    
    auto span8 = util::span(v.begin(), v.begin() + 8);
    std::cout << "Sum of 8 numbers: " << sum(span8) << std::endl;
    
    auto full_span = util::make_span(v);
    std::cout << "Sum of all numbers: " << sum(full_span) << std::endl;
    
  } // analyse()
  
}; // struct RangeDocumentation1TestClass


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(span_testcase) {
  
  using TestVector_t = std::vector<int>;
  
  TestVector_t ev;
  test_span<TestVector_t::iterator>(ev);
  test_const_span<TestVector_t::const_iterator>(ev);
  
  TestVector_t v3 { 1, 2, 3 };
  test_span<TestVector_t::iterator>(v3);
  test_const_span<TestVector_t::const_iterator>(v3);
  
  TestVector_t const cv4 { 1, 2, 3, 4 };
  test_span<TestVector_t::const_iterator>(cv4);
  test_const_span<TestVector_t::const_iterator>(cv4);
  
  
} // BOOST_AUTO_TEST_CASE(span_testcase)


//------------------------------------------------------------------------------
BOOST_AUTO_TEST_CASE(span_documentation_testcase) {
  
  // this only checks that the example compiles and runs
  RangeDocumentation1TestClass doc1;
  doc1.analyse();
  
} // BOOST_AUTO_TEST_CASE(span_documentation_testcase)


//------------------------------------------------------------------------------
