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
#include <memory> // std::unique_ptr<>
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

//
// C array container
//
static_assert(std::is_same_v<util::collection_value_t                <O      *>, O       >);
static_assert(std::is_same_v<util::collection_value_access_t         <O      *>, O      &>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<O      *>, O const&>);

static_assert(std::is_same_v<util::collection_value_t                <O const*>, O const >);
static_assert(std::is_same_v<util::collection_value_access_t         <O const*>, O const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<O const*>, O const&>);

static_assert(std::is_same_v<util::collection_value_t                <O[10U]  >, O       >);
static_assert(std::is_same_v<util::collection_value_access_t         <O[10U]  >, O      &>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<O[10U]  >, O const&>);

static_assert(std::is_same_v<util::collection_value_t                <std::unique_ptr<O      >>, O       >);
static_assert(std::is_same_v<util::collection_value_access_t         <std::unique_ptr<O      >>, O      &>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::unique_ptr<O      >>, O const&>);

static_assert(std::is_same_v<util::collection_value_t                <std::unique_ptr<O const>>, O const >);
static_assert(std::is_same_v<util::collection_value_access_t         <std::unique_ptr<O const>>, O const&>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::unique_ptr<O const>>, O const&>);

static_assert(std::is_same_v<util::collection_value_t                <std::unique_ptr<O[10U] >>, O       >);
static_assert(std::is_same_v<util::collection_value_access_t         <std::unique_ptr<O[10U] >>, O      &>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::unique_ptr<O[10U] >>, O const&>);

static_assert(std::is_same_v<util::collection_value_t                <std::unique_ptr<O[]    >>, O       >);
static_assert(std::is_same_v<util::collection_value_access_t         <std::unique_ptr<O[]    >>, O      &>);
static_assert(std::is_same_v<util::collection_value_constant_access_t<std::unique_ptr<O[]    >>, O const&>);

//
// util::collection_reference_t
//
static_assert(std::is_same_v<util::collection_reference_t<                       std::list<O> >, std::reference_wrapper<std::list<O>      >>);
static_assert(std::is_same_v<util::collection_reference_t<                       std::list<O>&>, std::reference_wrapper<std::list<O>      >>);
static_assert(std::is_same_v<util::collection_reference_t<std::reference_wrapper<std::list<O>>>, std::reference_wrapper<std::list<O>      >>);
static_assert(std::is_same_v<util::collection_reference_t<                       O*           >,                                  O*       >);
static_assert(std::is_same_v<util::collection_reference_t<                       O const*     >,                                  O const* >);
static_assert(std::is_same_v<util::collection_reference_t<                       O[10U]       >,                                  O*       >);


//------------------------------------------------------------------------------
void make_collection_reference_test() {
  
  /*
   * `std::unique_ptr<T[N]>` is a strange beast which really deals with pointers
   * to whole C arrays (`T(*)[N]`).
   * It should be understood that this is not the same as `T*`
   */
  constexpr std::size_t Size = 10U;
  
  std::list<O> list        {{ O{} }};
  O            array       [Size];
  O*           ptr         = array;
  O const*     cptr        = ptr;
  auto         refw        = std::ref(list);
  auto         crefw       = std::cref(list);
  auto const   uptr        = std::unique_ptr<O>(new O[Size]);
  auto         array_uptr  = std::unique_ptr<O[Size]>();
  auto         garray_uptr = std::unique_ptr<O[]>(new O[Size]);
  
  static_assert(std::is_same_v<decltype(util::make_collection_reference(list       )), decltype(std::ref(list))>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(array      )), decltype(ptr           )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(ptr        )), decltype(ptr           )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(cptr       )), decltype(cptr          )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(refw       )), decltype(refw          )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(crefw      )), decltype(crefw         )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(uptr       )), decltype(ptr           )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(array_uptr )), decltype(&array        )>);
  static_assert(std::is_same_v<decltype(util::make_collection_reference(garray_uptr)), decltype(ptr           )>);
  
  BOOST_CHECK_EQUAL(&(util::make_collection_reference(list       ).get().front()), &list.front()      );
  BOOST_CHECK_EQUAL(&(util::make_collection_reference(array      )[0]           ), &array[0]          );
  BOOST_CHECK_EQUAL(&(util::make_collection_reference(ptr        )[0]           ),  ptr               );
  BOOST_CHECK_EQUAL(&(util::make_collection_reference(cptr       )[0]           ),  cptr              );
  BOOST_CHECK_EQUAL(&(util::make_collection_reference(refw       ).get().front()), &list.front()      );
  BOOST_CHECK_EQUAL(&(util::make_collection_reference(crefw      ).get().front()), &list.front()      );
  BOOST_CHECK_EQUAL(&(util::make_collection_reference(uptr       )[0]           ),  uptr.get()        );
  BOOST_CHECK_EQUAL(&(util::make_collection_reference(array_uptr )[0]           ),  array_uptr.get()  );
  BOOST_CHECK_EQUAL(&(util::make_collection_reference(garray_uptr)[0]           ),  garray_uptr.get() );
  
} // make_collection_reference_test()


//------------------------------------------------------------------------------

BOOST_AUTO_TEST_CASE(NothingTestCase) {
  make_collection_reference_test();
} // BOOST_AUTO_TEST_CASE(NothingTestCase)

//------------------------------------------------------------------------------
