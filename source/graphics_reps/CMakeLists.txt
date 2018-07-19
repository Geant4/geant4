#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4graphics_reps
# Package: Geant4.src.G4graphics_reps
#
# CMakeLists.txt for single level library that may be build global or granular
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 108187 2018-01-18 16:03:33Z gcosmo $
#
#------------------------------------------------------------------------------

include(Geant4MacroLibraryTargets)

# Add allocation export symbol for the graphics_reps category
add_definitions(-DG4GREPS_ALLOC_EXPORT)

if(GEANT4_BUILD_GRANULAR_LIBS)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
else()
    GEANT4_GLOBAL_LIBRARY_TARGET(COMPONENTS sources.cmake)
endif()

