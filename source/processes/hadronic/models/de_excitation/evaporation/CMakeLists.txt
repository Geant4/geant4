#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4hadronic_deex_evaporation
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_deex.G4hadronic_deex_evaporation
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt,v 1.1 2010-09-29 18:59:18 bmorgan Exp $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

