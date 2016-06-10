#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hntools
# Package: Geant4.src.G4analysis.G4hntools
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 15/07/2013
#
# $Id: sources.cmake 92688 2015-09-14 07:01:13Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/g4tools/include)
include_directories(${CMAKE_SOURCE_DIR}/source/analysis/management/include)

#
# Define the Geant4 Module.
#

include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hntools
    HEADERS
        G4BaseHistoUtilities.hh
        G4MPIToolsManager.hh
        G4H1ToolsManager.hh
        G4H2ToolsManager.hh
        G4H3ToolsManager.hh
        G4P1ToolsManager.hh
        G4P2ToolsManager.hh
        G4ToolsAnalysisManager.hh
        G4ToolsAnalysisManager.icc
        G4ToolsAnalysisReader.hh
        G4ToolsAnalysisReader.icc
    SOURCES
        G4BaseHistoUtilities.cc
        G4H1ToolsManager.cc
        G4H2ToolsManager.cc
        G4H3ToolsManager.cc
        G4P1ToolsManager.cc
        G4P2ToolsManager.cc
        G4ToolsAnalysisManager.cc
        G4ToolsAnalysisReader.cc
   GRANULAR_DEPENDENCIES
        G4globman
        G4intercoms
        G4analysismng
    GLOBAL_DEPENDENCIES
        G4global
        G4intercoms
    LINK_LIBRARIES
)

# List any source specific properties here
