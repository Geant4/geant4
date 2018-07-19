#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4channeling
# Package: Geant4.src.G4processes.G4solidstate.G4channeling
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 06/04/2017
#
# $Id: CMakeLists.txt $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

