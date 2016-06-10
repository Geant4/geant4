#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4transportation
# Package: Geant4.src.G4processes.G4transportation
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 69887 2013-05-17 08:17:02Z gcosmo $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

