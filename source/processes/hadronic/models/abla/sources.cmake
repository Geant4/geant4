# - G4hadronic_abla module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_abla
  PUBLIC_HEADERS
    G4AblaDataDefs.hh
    G4AblaRandom.hh
    G4AblaInterface.hh
    G4AblaDataFile.hh
    G4Abla.hh
    G4AblaVirtualData.hh
  SOURCES
    G4AblaVirtualData.cc
    G4Abla.cc
    G4AblaRandom.cc
    G4AblaDataFile.cc
    G4AblaInterface.cc)

geant4_module_link_libraries(G4hadronic_abla
  PUBLIC
    G4globman
    G4hadronic_inclxx_utils
    G4hadronic_mgt
    G4hadronic_util
  PRIVATE
    G4hadronic_deex_handler
    G4heprandom
    G4ions
    G4partman)
