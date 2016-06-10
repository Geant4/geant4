#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4hadronic_proc
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_proc
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66499 2012-12-19 09:16:35Z gcosmo $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

