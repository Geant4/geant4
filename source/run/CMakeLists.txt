#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4run
# Package: Geant4.src.G4run
#
# CMakeLists.txt for single level library that may be build global or granular
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 68802 2013-04-05 13:54:41Z gcosmo $
#
#------------------------------------------------------------------------------

# Add allocation export symbol for the run category
add_definitions(-DG4RUN_ALLOC_EXPORT)

include(Geant4MacroLibraryTargets)
if(GEANT4_BUILD_GRANULAR_LIBS)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
else()
    GEANT4_GLOBAL_LIBRARY_TARGET(COMPONENTS sources.cmake)
endif()

