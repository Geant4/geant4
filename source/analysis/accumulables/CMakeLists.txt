#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4analysis
# Package: Geant4.src.G4analysis.G4analysismng
#
# CMakeLists.txt for building a single granular library.
#
# Created on : 15/07/2013
#
# $Id: CMakeLists.txt 99561 2016-09-27 07:04:12Z gcosmo $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

