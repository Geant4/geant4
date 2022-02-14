# G4biasing_mgt module build definition

# Define the Geant4 Module.
geant4_add_module(G4biasing_mgt
  PUBLIC_HEADERS
    G4VProcessPlacer.hh
    G4ProcessPlacer.hh
    G4BiasingAppliedCase.hh
    G4BiasingOperationManager.hh
    G4VBiasingInteractionLaw.hh
    G4VBiasingOperation.hh
    G4VBiasingOperator.hh
  SOURCES
    G4VProcessPlacer.cc
    G4ProcessPlacer.cc
    G4BiasingOperationManager.cc
    G4VBiasingOperation.cc
    G4VBiasingOperator.cc)

geant4_module_link_libraries(G4biasing_mgt
  PUBLIC
    G4globman
    G4track
  PRIVATE
    G4partman
    G4procman)
