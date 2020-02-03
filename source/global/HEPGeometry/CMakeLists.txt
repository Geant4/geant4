#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4hepgeometry
# Package: Geant4.src.G4global.G4hepgeometry
#
# CMakeLists.txt for building a header-only module.
#
# Generated on : 24/9/2010
#
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_HEADER_MODULE_TARGET(COMPONENT sources.cmake)
endif()

