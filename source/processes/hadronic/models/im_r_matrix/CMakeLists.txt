#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4had_im_r_matrix
# Package: Geant4.src.G4processes.G4hadronic.G4hadronic_models.G4had_im_r_matrix
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 69856 2013-05-16 13:32:36Z gcosmo $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

