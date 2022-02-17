cet_make(NO_DICTIONARY
         LIBRARIES
           ROOT::Core
           ROOT::Physics
         ) ### MIGRATE-ACTION-RECOMMENDED (migrate-3.05.03) - deprecated: use cet_make_library(), build_dictonary(), cet_plugin() with explicit source lists and plugin base types

art_dictionary(DICTIONARY_LIBRARIES larcorealg_GeoAlgo)


install_headers()
install_source()
