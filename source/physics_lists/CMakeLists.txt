#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4physicslists
# Package: Geant4.src.G4physicslists
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt,v 1.1 2010-09-29 18:50:46 bmorgan Exp $
#
#------------------------------------------------------------------------------

add_subdirectory(builders)
add_subdirectory(lists)

if(NOT GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4physicslists
        COMPONENTS
            builders/sources.cmake
            lists/sources.cmake
)
endif()

