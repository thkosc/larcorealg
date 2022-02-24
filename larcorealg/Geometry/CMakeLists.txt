cet_make_library(
  SOURCE AuxDetChannelMapAlg.cxx
         AuxDetGeo.cxx
         AuxDetGeometryCore.cxx
         AuxDetSensitiveGeo.cxx
         BoxBoundedGeo.cxx
         ChannelMapAlg.cxx
         ChannelMapStandardAlg.cxx
         CryostatGeo.cxx
         DriftPartitions.cxx
         GeometryBuilderStandard.cxx
         GeometryBuilderWireless.cxx
         GeometryCore.cxx
         GeoNodePath.cxx
         GeoObjectSorter.cxx
         GeoObjectSorterStandard.cxx
         GeoVectorLocalTransformation.cxx
         LocalTransformation.cxx
         OpDetGeo.cxx
         PlaneGeo.cxx
         StandaloneGeometrySetup.cxx
         TPCGeo.cxx
         WireGeo.cxx
  LIBRARIES
          messagefacility::MF_MessageLogger
          fhiclcpp::fhiclcpp
          cetlib::cetlib
          cetlib_except::cetlib_except
          ROOT::Core
          ROOT::Physics
          ROOT::Geom
          ROOT::GenVector
          CLHEP::Geometry
          PUBLIC
          larcoreobj::SimpleTypesAndConstants
        ) 

cet_make_library(LIBRARY_NAME geoVectorUtils INTERFACE
 SOURCE geo_vectors_fhicl.h geo_vectors_utils.h geo_vectors_utils_TVector.h
  LIBRARIES INTERFACE
          larcoreobj::headers
          fhiclcpp::fhiclcpp
          ROOT::GenVector
)

cet_make_library(LIBRARY_NAME GeometryIDmapper INTERFACE
 SOURCE  GeometryIDmapper.h
  LIBRARIES INTERFACE
          larcoreobj::SimpleTypesAndConstants
)

install_headers(SUBDIRS "details")
install_fhicl(SUBDIRS "details")
install_source(SUBDIRS "details")

add_subdirectory(gdml)
