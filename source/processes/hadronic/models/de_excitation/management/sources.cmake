# - G4hadronic_deex_management module build definition

# Define the Geant4 Module.
geant4_add_module(G4hadronic_deex_management
  PUBLIC_HEADERS
    G4DeexParametersMessenger.hh
    G4DeexPrecoParameters.hh
    G4LevelManager.hh
    G4LevelReader.hh
    G4NuclearLevelData.hh
    G4NucLevel.hh
    G4VEmissionProbability.hh
    G4VEvaporationChannel.hh
    G4VEvaporationFactory.hh
  SOURCES
    G4DeexParametersMessenger.cc
    G4DeexPrecoParameters.cc
    G4LevelManager.cc
    G4LevelReader.cc
    G4NuclearLevelData.cc
    G4NucLevel.cc
    G4VEmissionProbability.cc
    G4VEvaporationChannel.cc
    G4VEvaporationFactory.cc)

geant4_module_link_libraries(G4hadronic_deex_management
  PUBLIC
    G4globman
    G4hadronic_util
    G4intercoms
  PRIVATE
    G4hadronic_deex_util
    G4heprandom
    G4materials
    G4partman)
