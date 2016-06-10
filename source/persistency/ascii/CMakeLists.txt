#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4geomtext
# Package: Geant4.src.G4persistency.G4geomtext
#
# CMakeLists.txt for building a single granular library.
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66363 2012-12-18 09:12:54Z gcosmo $
#
#------------------------------------------------------------------------------

if(GEANT4_BUILD_GRANULAR_LIBS)
    include(Geant4MacroLibraryTargets)
    GEANT4_GRANULAR_LIBRARY_TARGET(COMPONENT sources.cmake)
endif()

