#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4geomdivision
# Package: Geant4.src.G4geometry.G4geomdivision
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 69784 2013-05-15 09:16:06Z gcosmo $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

