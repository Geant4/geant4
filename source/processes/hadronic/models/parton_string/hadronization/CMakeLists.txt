#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4had_string_frag
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_string.G4had_string_frag
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66892 2013-01-17 10:57:59Z gunter $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

