cet_test(NameSelector_test USE_BOOST_UNIT
  SOURCE NameSelector_test.cxx
  LIBRARIES PRIVATE
  larcorealg::TestUtils
  canvas::canvas
)

cet_test(ProviderList_test USE_BOOST_UNIT
  LIBRARIES PRIVATE
  larcorealg::ProviderTestHelpers
)

cet_test(StopWatch_test)
