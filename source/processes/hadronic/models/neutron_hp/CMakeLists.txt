#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4had_neu_hp
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_neu_hp
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 86780 2014-11-18 08:38:07Z gcosmo $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

