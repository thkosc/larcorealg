cet_make_library(SOURCE
  makeValueIndex.h
  DebugUtils.cxx
  SortByPointers.h
  LIBRARIES PUBLIC
  cetlib_except::cetlib_except
)

cet_make_library(LIBRARY_NAME Math_Functor INTERFACE
  SOURCE quiet_Math_Functor.h
  LIBRARIES INTERFACE
  ROOT::MathCore
)

cet_make_library(LIBRARY_NAME ParticleFilters INTERFACE
  SOURCE ParticleFilters.h
  LIBRARIES INTERFACE
  ROOT::Geom
  ROOT::Matrix
  ROOT::Physics
)

install_headers()
install_source()
