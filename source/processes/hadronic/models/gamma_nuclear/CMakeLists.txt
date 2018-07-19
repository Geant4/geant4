#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4had_gamm_nuclear
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4lepto_nuclear
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 09/Oct/2017
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

