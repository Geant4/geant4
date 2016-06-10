#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4gflash
# Package: Geant4.src.G4parmodels.G4gflash
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 69796 2013-05-15 13:26:12Z gcosmo $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

