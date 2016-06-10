#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4analysis
# Package: Geant4.src.G4analysis.G4root
#
# CMakeLists.txt for building a single granular library.
#
# Created on : 15/07/2013
#
# $Id: CMakeLists.txt 72292 2013-07-15 12:01:43Z ihrivnac $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

