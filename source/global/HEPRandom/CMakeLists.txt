#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : 
# Package: Geant4.src.G4global.
#
# CMakeLists.txt for building a header-only module.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt,v 1.1 2010-09-29 18:44:25 bmorgan Exp $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_HEADER_MODULE_TARGET(COMPONENT sources.cmake)
endif()

