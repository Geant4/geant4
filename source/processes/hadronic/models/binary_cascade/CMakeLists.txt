#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4hadronic_binary
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_binary
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66803 2013-01-12 15:49:41Z gcosmo $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

