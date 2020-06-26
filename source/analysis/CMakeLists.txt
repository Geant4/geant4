#------------------------------------------------------------------------------
# Module : G4geometry
# Package: Geant4.src.G4analysis
#------------------------------------------------------------------------------

if(GEANT4_USE_FREETYPE)
  add_definitions(-DTOOLS_USE_FREETYPE)
  add_definitions(-DTOOLS_USE_EXPAT)
endif()

# g4tools needs additional CMake setup
add_subdirectory(g4tools)

set(G4analysis_COMPONENTS
  accumulables/sources.cmake
  csv/sources.cmake
  factory/sources.cmake
  g4tools/sources.cmake
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
  COMPONENTS ${G4analysis_COMPONENTS}
  )

#-----------------------------------------------------------------------
# Install fonts
#
if(GEANT4_USE_FREETYPE)
  install(DIRECTORY fonts/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/fonts
    PATTERN ".svn" EXCLUDE
  )
endif()
