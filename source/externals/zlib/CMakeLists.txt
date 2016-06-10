#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4expat
# Package: Geant4.src.G4externals.G4zlib
#
# CMakeLists.txt for single level library that may be build global or granular
#
# Created on: 03/05/2013
#
# $Date$
# $Revision$
# $Author$
#
#------------------------------------------------------------------------------

include(Geant4MacroLibraryTargets)
if(GEANT4_BUILD_GRANULAR_LIBS)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
else()
    GEANT4_GLOBAL_LIBRARY_TARGET(COMPONENTS sources.cmake)
endif()
