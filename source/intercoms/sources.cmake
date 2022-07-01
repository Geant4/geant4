# - G4intercoms module build definition

geant4_add_module(G4intercoms
  PUBLIC_HEADERS
    G4AnyMethod.hh
    G4AnyType.hh
    G4GenericMessenger.hh
    G4LocalThreadCoutMessenger.hh
    G4ProfilerMessenger.hh
    G4UIaliasList.hh
    G4UIbatch.hh
    G4UIbridge.hh
    G4UIcmdWith3Vector.hh
    G4UIcmdWith3VectorAndUnit.hh
    G4UIcmdWithABool.hh
    G4UIcmdWithADouble.hh
    G4UIcmdWithADoubleAndUnit.hh
    G4UIcmdWithAString.hh
    G4UIcmdWithAnInteger.hh
    G4UIcmdWithALongInt.hh
    G4UIcmdWithoutParameter.hh
    G4UIcommand.hh
    G4UIcommandStatus.hh
    G4UIcommandTree.hh
    G4UIcontrolMessenger.hh
    G4UIdirectory.hh
    G4UImanager.hh
    G4UImessenger.hh
    G4UIparameter.hh
    G4UIsession.hh
    G4UItokenNum.hh
    G4UnitsMessenger.hh
    G4VFlavoredParallelWorld.hh
    G4VGlobalFastSimulationManager.hh
    icomsdefs.hh
  SOURCES
    G4LocalThreadCoutMessenger.cc
    G4GenericMessenger.cc
    G4ProfilerMessenger.cc
    G4UIaliasList.cc
    G4UIbatch.cc
    G4UIbridge.cc
    G4UIcmdWith3Vector.cc
    G4UIcmdWith3VectorAndUnit.cc
    G4UIcmdWithABool.cc
    G4UIcmdWithADouble.cc
    G4UIcmdWithADoubleAndUnit.cc
    G4UIcmdWithAString.cc
    G4UIcmdWithAnInteger.cc
    G4UIcmdWithALongInt.cc
    G4UIcmdWithoutParameter.cc
    G4UIcommand.cc
    G4UIcommandTree.cc
    G4UIcontrolMessenger.cc
    G4UIdirectory.cc
    G4UImanager.cc
    G4UImessenger.cc
    G4UIparameter.cc
    G4UIsession.cc
    G4UnitsMessenger.cc
    G4VGlobalFastSimulationManager.cc)

geant4_module_compile_definitions(G4intercoms PRIVATE G4ICOMS_ALLOC_EXPORT)

geant4_module_link_libraries(G4intercoms PUBLIC G4globman PRIVATE ${timemory_LIBRARIES})
