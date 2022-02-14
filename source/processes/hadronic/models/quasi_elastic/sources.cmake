# - G4hadronic_quasi_elastic module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_quasi_elastic
  PUBLIC_HEADERS G4QuasiElRatios.hh
  SOURCES G4QuasiElRatios.cc)

geant4_module_link_libraries(G4hadronic_quasi_elastic
  PUBLIC
    G4globman
    G4hadronic_xsect
    G4hepgeometry
    G4heprandom
  PRIVATE
    G4baryons
    G4ions)
