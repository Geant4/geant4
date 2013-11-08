#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4had_fission
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_fission
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 08/11/2013
#
# $Id:$
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

