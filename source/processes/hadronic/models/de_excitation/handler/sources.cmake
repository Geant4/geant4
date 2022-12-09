# - G4hadronic_deex_handler module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_deex_handler
  PUBLIC_HEADERS G4ExcitationHandler.hh
  SOURCES G4ExcitationHandler.cc)

geant4_module_link_libraries(G4hadronic_deex_handler
  PUBLIC
    G4globman
    G4hadronic_deex_management
    G4hadronic_util
    G4materials
    G4partman
  PRIVATE
    G4hadronic_deex_evaporation
    G4hadronic_deex_fermi_breakup
    G4hadronic_deex_multifragmentation
    G4hadronic_deex_photon_evaporation
    G4hepgeometry
    G4leptons
    G4baryons
    G4procman)
