#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4scoring
# Package: Geant4.src.G4processes.G4scoring
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 69885 2013-05-17 07:49:30Z gcosmo $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

