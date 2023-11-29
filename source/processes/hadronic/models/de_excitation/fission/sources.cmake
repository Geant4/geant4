# - G4hadronic_deex_fission module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_deex_fission
  PUBLIC_HEADERS
    G4CompetitiveFission.hh
    G4EvaporationLevelDensityParameter.hh
    G4FissionBarrier.hh
    G4FissionLevelDensityParameter.hh
    G4FissionLevelDensityParameterINCLXX.hh
    G4FissionParameters.hh
    G4FissionProbability.hh
    G4VFissionBarrier.hh
  SOURCES
    G4CompetitiveFission.cc
    G4EvaporationLevelDensityParameter.cc
    G4FissionBarrier.cc
    G4FissionLevelDensityParameter.cc
    G4FissionLevelDensityParameterINCLXX.cc
    G4FissionParameters.cc
    G4FissionProbability.cc
    G4VFissionBarrier.cc)

geant4_module_link_libraries(G4hadronic_deex_fission
  PUBLIC
    G4globman
    G4hadronic_deex_management
    G4hadronic_deex_util
    G4hadronic_util
  PRIVATE
    G4heprandom
    G4partman)
