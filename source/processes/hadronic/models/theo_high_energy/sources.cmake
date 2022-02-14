# - G4had_theo_max module build definition

# Define the Geant4 Module.
geant4_add_module(G4had_theo_max
  PUBLIC_HEADERS
    G4CRCoalescence.hh
    G4QuasiElasticChannel.hh
    G4TheoFSGenerator.hh
  SOURCES
    G4CRCoalescence.cc
    G4QuasiElasticChannel.cc
    G4TheoFSGenerator.cc)

geant4_module_link_libraries(G4had_theo_max
  PUBLIC
    G4globman
    G4hadronic_mgt
    G4hadronic_util
  PRIVATE
    G4baryons
    G4hadronic_quasi_elastic
    G4hepgeometry
    G4ions
    G4partman)
