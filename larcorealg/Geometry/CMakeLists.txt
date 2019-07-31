cet_make(LIBRARIES
          ${MF_MESSAGELOGGER}
          ${FHICLCPP}
          cetlib
          cetlib_except
          ${ROOT_CORE}
          ${ROOT_PHYSICS}
          ${ROOT_GEOM}
          ${ROOT_GENVECTOR}
        )

install_headers()
install_fhicl()
install_source()

add_subdirectory(gdml)
