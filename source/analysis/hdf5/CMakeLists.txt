#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4analysis
# Package: Geant4.src.G4analysis.G4xml
#
# CMakeLists.txt for building a single granular library.
#
# Created on : 27/07/2017
#
# $Id: CMakeLists.txt 105426 2017-07-27 09:10:08Z ihrivnac $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

