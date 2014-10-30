#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4geomUsolids
# Package: Geant4.src.G4externals.G4geomUSolids
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66356 2012-12-18 09:02:32Z gcosmo $
#
#------------------------------------------------------------------------------

include(Geant4MacroLibraryTargets)
if(GEANT4_BUILD_GRANULAR_LIBS)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
else()
    GEANT4_GLOBAL_LIBRARY_TARGET(COMPONENTS sources.cmake)
endif()
