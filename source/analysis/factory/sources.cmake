# - G4analysisfac module build definition

# Optional additional links
if(GEANT4_USE_HDF5)
  set(G4analysisfac_G4hdf5 G4hdf5)
endif()

# Define the Geant4 Module.
geant4_add_module(G4analysisfac
  PUBLIC_HEADERS
    G4AnalysisManager.hh
    G4GenericAnalysisManager.hh
    G4GenericAnalysisManager.icc
    G4GenericFileManager.hh
    G4GenericFileManager.icc
  SOURCES
    G4GenericAnalysisManager.cc
    G4GenericFileManager.cc)

geant4_module_link_libraries(G4analysisfac
  PUBLIC G4analysismng G4hntools G4globman
  PRIVATE G4csv G4root G4xml ${G4analysisfac_G4hdf5})
