#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4magneticfield
# Package: Geant4.src.G4geometry.G4magneticfield
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66356 2012-12-18 09:02:32Z gcosmo $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

#if(GEANT4_ENABLE_TESTING)
#  add_subdirectory(test)
#endif()
