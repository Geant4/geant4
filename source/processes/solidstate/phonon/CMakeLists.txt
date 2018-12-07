#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4phonon
# Package: Geant4.src.G4processes.G4solidstate.G4phonon
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/10/2013
#
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

