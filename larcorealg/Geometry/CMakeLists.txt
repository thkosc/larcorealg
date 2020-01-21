cet_make(
  SUBDIRS
    "details"
  LIBRARIES
          ${MF_MESSAGELOGGER}
          ${FHICLCPP}
          cetlib
          cetlib_except
          ${ROOT_CORE}
          ${ROOT_PHYSICS}
          ${ROOT_GEOM}
          ${ROOT_GENVECTOR}
        )

install_headers(SUBDIRS "details")
install_fhicl(SUBDIRS "details")
install_source(SUBDIRS "details")

add_subdirectory(gdml)
