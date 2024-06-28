# - G4hadronic_nudex module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_nudex
  PUBLIC_HEADERS
    G4NuDEXNeutronCaptureModel.hh
  PRIVATE_HEADERS
    G4NuDEXInternalConversion.hh
    G4NuDEXLevelDensity.hh
    G4NuDEXPSF.hh
    G4NuDEXRandom.hh
    G4NuDEXStatisticalNucleus.hh
  SOURCES
    G4NuDEXInternalConversion.cc
    G4NuDEXLevelDensity.cc
    G4NuDEXNeutronCaptureModel.cc
    G4NuDEXPSF.cc
    G4NuDEXRandom.cc
    G4NuDEXStatisticalNucleus.cc)

geant4_module_link_libraries(G4hadronic_nudex
  PUBLIC
    G4globman
    G4hadronic_mgt
    G4hadronic_util
  PRIVATE
    G4bosons
    G4ions
    G4leptons
    G4partman
    G4hadronic_deex_management
    G4hadronic_deex_photon_evaporation
    G4hepgeometry
    G4heprandom)
