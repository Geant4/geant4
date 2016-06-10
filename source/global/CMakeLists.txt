#-----------------------------------------------------------------------
# CMakeLists.txt
# Module : G4global
# Package: Geant4.src.G4global
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 69106 2013-04-18 13:05:15Z gcosmo $
#
#-----------------------------------------------------------------------

# Add allocation export symbol for the global category
add_definitions(-DG4GLOB_ALLOC_EXPORT)

add_subdirectory(HEPGeometry)
add_subdirectory(HEPNumerics)
add_subdirectory(HEPRandom)
add_subdirectory(management)

# Geant4 global library
if(NOT GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4global
        COMPONENTS
            HEPGeometry/sources.cmake
            HEPNumerics/sources.cmake
            HEPRandom/sources.cmake
            management/sources.cmake
)
endif()

