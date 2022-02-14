# - G4hadronic_deex_photon_evaporation module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_deex_photon_evaporation
  PUBLIC_HEADERS
    G4GammaTransition.hh
    G4NeutronRadCapture.hh
    G4PhotonEvaporation.hh
    G4PolarizationTransition.hh
    G4VGammaTransition.hh
  SOURCES
    G4GammaTransition.cc
    G4NeutronRadCapture.cc
    G4PhotonEvaporation.cc
    G4PolarizationTransition.cc)

geant4_module_link_libraries(G4hadronic_deex_photon_evaporation
  PUBLIC
    G4globman
    G4hadronic_deex_management
    G4hadronic_mgt
    G4hadronic_util
    G4hepgeometry
  PRIVATE
    G4bosons
    G4heprandom
    G4ions
    G4leptons
    G4materials
    G4partman)
