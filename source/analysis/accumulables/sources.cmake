# - G4accumulables module build definition

# Define the Geant4 Module.
geant4_add_module(G4accumulables
  PUBLIC_HEADERS
    G4MergeMode.hh
    G4Accumulable.hh
    G4AccumulableManager.hh
    G4AccumulableManager.icc
    G4AccValue.hh
    G4AccValue.icc
    G4AccArray.hh
    G4AccArray.icc
    G4AccMap.hh
    G4AccMap.icc
    G4AccUnorderedMap.hh
    G4AccUnorderedMap.icc
    G4AccVector.hh
    G4AccVector.icc
    G4AccType.hh
    G4VAccumulable.hh
    G4VAccumulable.icc
  SOURCES
    G4MergeMode.cc
    G4AccumulableManager.cc)

geant4_module_link_libraries(G4accumulables PUBLIC G4globman)
