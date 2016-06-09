#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4interfaces
# Package: Geant4.src.G4interfaces
#
# Top level CMakeLists.txt for a Geant4 Global Library
# or collection of sublibraries
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt,v 1.1 2010-09-29 18:45:30 bmorgan Exp $
#
#------------------------------------------------------------------------------

add_subdirectory(GAG)
add_subdirectory(basic)
add_subdirectory(common)

if(NOT GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GLOBAL_LIBRARY_TARGET(NAME G4interfaces
        COMPONENTS
            GAG/sources.cmake
            basic/sources.cmake
            common/sources.cmake
)
endif()

