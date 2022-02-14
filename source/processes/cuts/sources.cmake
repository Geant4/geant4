# G4cuts module build definition

# Define the Geant4 Module.
geant4_add_module(G4cuts
  PUBLIC_HEADERS
    G4MCCIndexConversionTable.hh
    G4MaterialCutsCouple.hh
    G4PhysicsTableHelper.hh
    G4ProductionCuts.hh
    G4ProductionCutsTable.hh
    G4ProductionCutsTableMessenger.hh
    G4RToEConvForElectron.hh
    G4RToEConvForGamma.hh
    G4RToEConvForPositron.hh
    G4RToEConvForProton.hh
    G4VRangeToEnergyConverter.hh
  SOURCES
    G4MCCIndexConversionTable.cc
    G4MaterialCutsCouple.cc
    G4PhysicsTableHelper.cc
    G4ProductionCuts.cc
    G4ProductionCutsTable.cc
    G4ProductionCutsTableMessenger.cc
    G4RToEConvForElectron.cc
    G4RToEConvForGamma.cc
    G4RToEConvForPositron.cc
    G4RToEConvForProton.cc
    G4VRangeToEnergyConverter.cc)

geant4_module_link_libraries(G4cuts
  PUBLIC
    G4geometrymng
    G4globman
    G4intercoms
    G4materials
    G4partman)
