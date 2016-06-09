#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4hadronic_coherent_elastic
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4hadronic_coherent_elastic
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt,v 1.1 2010-09-29 18:58:37 bmorgan Exp $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

