#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4OpenGL
# Package: Geant4.src.G4visualization.G4OpenGL
#
# CMakeLists.txt for single level library that may be build global or granular
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt,v 1.1 2010-09-29 19:12:06 bmorgan Exp $
#
#------------------------------------------------------------------------------

include(Geant4MacroLibraryTargets)
if(GEANT4_BUILD_GRANULAR_LIBS)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
else()
    GEANT4_GLOBAL_LIBRARY_TARGET(COMPONENTS sources.cmake)
endif()

