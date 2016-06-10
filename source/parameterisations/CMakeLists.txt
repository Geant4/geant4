#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4parmodels
# Package: Geant4.src.G4parmodels
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66892 2013-01-17 10:57:59Z gunter $
#
#------------------------------------------------------------------------------

add_subdirectory(gflash)

if(NOT GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4parmodels
        COMPONENTS
            gflash/sources.cmake
)
endif()

