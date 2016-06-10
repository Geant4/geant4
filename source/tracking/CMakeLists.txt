#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4tracking
# Package: Geant4.src.G4tracking
#
# CMakeLists.txt for single level library that may be build global or granular
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66821 2013-01-12 16:22:43Z gcosmo $
#
#------------------------------------------------------------------------------

# Add allocation export symbol for the tracking category
add_definitions(-DG4TRACKING_ALLOC_EXPORT)

include(Geant4MacroLibraryTargets)
if(GEANT4_BUILD_GRANULAR_LIBS)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
else()
    GEANT4_GLOBAL_LIBRARY_TARGET(COMPONENTS sources.cmake)
endif()

