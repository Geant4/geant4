# - G4hadronic_em_dissociation module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_em_dissociation
  PUBLIC_HEADERS G4EMDissociation.hh
  SOURCES G4EMDissociation.cc)

geant4_module_link_libraries(G4hadronic_em_dissociation
  PUBLIC
    G4globman
    G4hadronic_deex_handler
    G4hadronic_mgt
    G4hadronic_util
    G4hadronic_xsect
  PRIVATE
    G4baryons
    G4hepgeometry
    G4heprandom
    G4partman)
