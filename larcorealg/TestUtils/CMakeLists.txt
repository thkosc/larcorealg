cet_make_library(
         SOURCE NameSelector.cxx
         LIBRARIES canvas::canvas
                   cetlib_except::cetlib_except)

install_headers()
install_source()
