#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4empii
# Package: Geant4.src.G4processes.G4electromagnetic.G4empii
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 19/11/2010
#
# $Id: CMakeLists.txt,v 1.1 2010-11-19 16:05:03 gcosmo Exp $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

