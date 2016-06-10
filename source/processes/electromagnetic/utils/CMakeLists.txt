#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4emutils
# Package: Geant4.src.G4processes.G4electromagnetic.G4emutils
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66791 2013-01-12 15:24:24Z gcosmo $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

