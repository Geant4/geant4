# G4biasing_imp module build definition

# Define the Geant4 Module.
geant4_add_module(G4biasing_imp
  PUBLIC_HEADERS
    G4GeometrySampler.hh
    G4ImportanceConfigurator.hh
    G4ImportanceProcess.hh
    G4PlaceOfAction.hh
    G4SamplingPostStepAction.hh
    G4VSampler.hh
    G4VSamplerConfigurator.hh
    G4WeightCutOffConfigurator.hh
    G4WeightCutOffProcess.hh
    G4WeightWindowConfigurator.hh
    G4WeightWindowProcess.hh
  SOURCES
    G4GeometrySampler.cc
    G4ImportanceConfigurator.cc
    G4ImportanceProcess.cc
    G4SamplingPostStepAction.cc
    G4VSampler.cc
    G4VSamplerConfigurator.cc
    G4WeightCutOffConfigurator.cc
    G4WeightCutOffProcess.cc
    G4WeightWindowConfigurator.cc
    G4WeightWindowProcess.cc)

geant4_module_link_libraries(G4biasing_imp
  PUBLIC
    G4biasing_mgt
    G4geombias
    G4geometrymng
    G4globman
    G4magneticfield
    G4navigation
    G4procman
    G4transportation
  PRIVATE
    G4track)
