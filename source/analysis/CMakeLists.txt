# - G4analysis category build

if(GEANT4_USE_FREETYPE)
  add_definitions(-DTOOLS_USE_FREETYPE)
endif()

set(G4analysis_COMPONENTS
  accumulables/sources.cmake
  csv/sources.cmake
  factory/sources.cmake
  hntools/sources.cmake
  management/sources.cmake
  root/sources.cmake
  xml/sources.cmake
)

if(GEANT4_USE_HDF5)
  add_definitions(-DTOOLS_USE_HDF5)
  list(APPEND G4analysis_COMPONENTS hdf5/sources.cmake)
endif()


geant4_global_library_target(NAME G4analysis
  COMPONENTS ${G4analysis_COMPONENTS})
