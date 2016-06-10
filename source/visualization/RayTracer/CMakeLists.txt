#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4RayTracer
# Package: Geant4.src.G4visualization.G4RayTracer
#
# CMakeLists.txt for single level library that may be build global or granular
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66373 2012-12-18 09:41:34Z gcosmo $
#
#------------------------------------------------------------------------------

# Add allocation export symbol for the RayTracer module
# NB: we don't add this at the category level because RayTracer is
# the only vis driver that uses exports at the moment. We can move the
# add_definitions higher later if needed.
add_definitions(-DG4VIS_ALLOC_EXPORT)


include(Geant4MacroLibraryTargets)
if(GEANT4_BUILD_GRANULAR_LIBS)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
else()
    GEANT4_GLOBAL_LIBRARY_TARGET(COMPONENTS sources.cmake)
endif()

