#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4geometry
# Package: Geant4.src.G4analysis
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 105747 2017-08-16 13:16:13Z gcosmo $
#
#------------------------------------------------------------------------------

if (GEANT4_USE_FREETYPE)
  add_definitions(-DTOOLS_USE_FREETYPE)
  add_definitions(-DTOOLS_USE_EXPAT)
endif()

if (GEANT4_USE_HDF5)
  add_definitions(-DTOOLS_USE_HDF5)
endif()

add_subdirectory(accumulables)
add_subdirectory(csv)
add_subdirectory(g4tools)
add_subdirectory(hntools)
add_subdirectory(management)
add_subdirectory(root)
add_subdirectory(xml)
if (GEANT4_USE_HDF5)
  add_subdirectory(hdf5)
endif()

if(NOT GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    if (GEANT4_USE_HDF5)
      GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4analysis
          COMPONENTS
              accumulables/sources.cmake
              csv/sources.cmake
              g4tools/sources.cmake
              hdf5/sources.cmake
              hntools/sources.cmake
              management/sources.cmake
              root/sources.cmake
              xml/sources.cmake
    )
    else()
      GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4analysis
          COMPONENTS
              accumulables/sources.cmake
              csv/sources.cmake
              g4tools/sources.cmake
              hntools/sources.cmake
              management/sources.cmake
              root/sources.cmake
              xml/sources.cmake
      )
    endif()
endif()

#-----------------------------------------------------------------------
# Install fonts
#
if (GEANT4_USE_FREETYPE)
  install(DIRECTORY fonts/
    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Geant4-${Geant4_VERSION}/fonts
    PATTERN ".svn" EXCLUDE
  )
endif()



