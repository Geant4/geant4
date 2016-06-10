#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4physlist_ctor_factory
# Package: Geant4.src.G4physicslists.G4phys_ctors.G4phys_ctor_em
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 12-Jan-2013
#
# $Id: CMakeLists.txt 66704 2013-01-10 18:20:17Z gunter $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()


