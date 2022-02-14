# - G4hadronic_abrasion module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_abrasion
  PUBLIC_HEADERS
    G4NuclearAbrasionGeometry.hh
    G4WilsonAbrasionModel.hh
  SOURCES
    G4NuclearAbrasionGeometry.cc
    G4WilsonAbrasionModel.cc)

geant4_module_link_libraries(G4hadronic_abrasion
  PUBLIC
    G4globman
    G4hadronic_ablation
    G4hadronic_deex_handler
    G4hadronic_mgt
    G4hadronic_util
    G4track
  PRIVATE
    G4hadronic_deex_evaporation
    G4hepgeometry
    G4heprandom
    G4partman)
