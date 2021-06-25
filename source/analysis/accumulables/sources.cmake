# - G4accumulables module build definition

# Define the Geant4 Module.
geant4_add_module(G4accumulables
  PUBLIC_HEADERS
    G4MergeMode.hh
    G4AccumulableManager.hh
    G4AccumulableManager.icc
    G4Accumulable.hh
    G4Accumulable.icc
    G4VAccumulable.hh
    G4VAccumulable.icc
  SOURCES
    G4MergeMode.cc
    G4AccumulableManager.cc)

geant4_module_link_libraries(G4accumulables PUBLIC G4globman)
