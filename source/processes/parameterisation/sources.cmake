# - G4parameterisation module build definition

# Define the Geant4 Module.
geant4_add_module(G4parameterisation
  PUBLIC_HEADERS
    G4FastHit.hh
    G4FastSimHitMaker.hh
    G4FastSimulationHelper.hh
    G4FastSimulationManager.hh
    G4FastSimulationManagerProcess.hh
    G4FastSimulationProcessType.hh
    G4FastSimulationVector.hh
    G4FastSimulationVector.icc
    G4FastStep.hh
    G4FastStep.icc
    G4FastTrack.hh
    G4GlobalFastSimulationManager.hh
    G4VFastSimSensitiveDetector.hh
    G4VFastSimulationModel.hh
  PRIVATE_HEADERS
    G4FastSimulationMessenger.hh
  SOURCES
    G4FastSimHitMaker.cc
    G4FastSimulationHelper.cc
    G4FastSimulationManager.cc
    G4FastSimulationManagerProcess.cc
    G4FastSimulationMessenger.cc
    G4FastStep.cc
    G4FastTrack.cc
    G4GlobalFastSimulationManager.cc
    G4VFastSimulationModel.cc)

geant4_module_link_libraries(G4parameterisation
  PUBLIC
    G4detector
    G4geometrymng
    G4globman
    G4hepgeometry
    G4intercoms
    G4magneticfield
    G4navigation
    G4partman
    G4procman
    G4track
  PRIVATE
    G4materials
    G4volumes)
