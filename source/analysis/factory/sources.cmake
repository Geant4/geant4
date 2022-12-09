# - G4analysisfac module build definition

# Define the Geant4 Module.
geant4_add_module(G4analysisfac
  PUBLIC_HEADERS
    G4AnalysisManager.hh
    G4GenericAnalysisManager.hh
    G4GenericAnalysisManager.icc
    G4GenericAnalysisMessenger.hh
    G4GenericFileManager.hh
    G4GenericFileManager.icc
  SOURCES
    G4GenericAnalysisManager.cc
    G4GenericAnalysisMessenger.cc
    G4GenericFileManager.cc)

geant4_module_link_libraries(G4analysisfac
  PUBLIC G4analysismng G4hntools G4globman G4intercoms
  PRIVATE G4csv G4root G4xml)

# HDF5, if enabled
if(GEANT4_USE_HDF5)
  geant4_module_compile_definitions(G4analysisfac PUBLIC TOOLS_USE_HDF5)
  geant4_module_link_libraries(G4analysisfac PRIVATE G4hdf5)
endif()
