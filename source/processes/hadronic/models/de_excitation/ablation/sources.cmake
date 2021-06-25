# - G4hadronic_ablation module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_ablation
  PUBLIC_HEADERS G4WilsonAblationModel.hh
  SOURCES G4WilsonAblationModel.cc)

geant4_module_link_libraries(G4hadronic_ablation
  PUBLIC
    G4globman
    G4hadronic_deex_evaporation
    G4hadronic_util
    G4partman
  PRIVATE
    G4baryons
    G4hadronic_deex_management
    G4hadronic_deex_photon_evaporation
    G4hepgeometry
    G4heprandom
    G4ions)
