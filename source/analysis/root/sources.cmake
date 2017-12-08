#------------------------------------------------------------------------------
# sources.cmake
# Module : G4csv
# Package: Geant4.src.G4analysis.G4root
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 15/07/2013
#
# $Id: sources.cmake 105338 2017-07-21 09:14:27Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${ZLIB_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/g4tools/include)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/hntools/include)
include_directories(${CMAKE_SOURCE_DIR}/source/run/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4root
    HEADERS
        G4RootAnalysisManager.hh
        G4RootAnalysisManager.icc
        G4RootAnalysisReader.hh
        G4RootAnalysisReader.icc
        G4RootFileManager.hh
        G4RootMainNtupleManager.hh
        G4RootNtupleManager.hh
        G4RootPNtupleDescription.hh
        G4RootPNtupleManager.hh
        G4RootRFileManager.hh
        G4RootRNtupleManager.hh
       g4root_defs.hh
        g4root.hh
    SOURCES
        G4RootAnalysisManager.cc
        G4RootAnalysisReader.cc
        G4RootFileManager.cc
        G4RootMainNtupleManager.cc
        G4RootNtupleManager.cc
        G4RootPNtupleManager.cc
        G4RootRFileManager.cc
        G4RootRNtupleManager.cc
    GRANULAR_DEPENDENCIES
        G4globman
        G4intercoms
        G4analysismng
        G4hntools
    GLOBAL_DEPENDENCIES
        G4global
        G4intercoms
    LINK_LIBRARIES
        ${ZLIB_LIBRARIES}
)

# List any source specific properties here
