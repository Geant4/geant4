#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4specsolids
# Package: Geant4.src.G4geometry..G4specsolids
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 67011 2013-01-29 16:17:41Z gcosmo $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

