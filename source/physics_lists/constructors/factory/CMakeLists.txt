#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4physlist_ctor_factory
# Package: Geant4.src.G4physicslists.G4phys_ctors.G4phys_ctor_factory
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66241 2012-12-13 18:34:42Z gunter $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()


