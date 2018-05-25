#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4event
# Package: Geant4.src.G4event
#
# CMakeLists.txt for single level library that may be build global or granular
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 110058 2018-05-15 09:02:41Z gcosmo $
#
#------------------------------------------------------------------------------

# Add allocation export symbol for the event category
add_definitions(-DG4EVENT_ALLOC_EXPORT)
########add_definitions(-DG4_USESMARTSTACK)

include(Geant4MacroLibraryTargets)
if(GEANT4_BUILD_GRANULAR_LIBS)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
else()
    GEANT4_GLOBAL_LIBRARY_TARGET(COMPONENTS sources.cmake)
endif()

