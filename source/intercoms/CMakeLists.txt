#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4intercoms
# Package: Geant4.src.G4intercoms
#
# CMakeLists.txt for single level library that may be build global or granular
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 106172 2017-09-15 13:03:57Z gcosmo $
#
#------------------------------------------------------------------------------

include(Geant4MacroLibraryTargets)

# Add allocation export symbol for the geometry category
add_definitions(-DG4ICOMS_ALLOC_EXPORT)

if(GEANT4_BUILD_GRANULAR_LIBS)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
else()
    GEANT4_GLOBAL_LIBRARY_TARGET(COMPONENTS sources.cmake)
endif()

