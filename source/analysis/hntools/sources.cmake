#------------------------------------------------------------------------------
# Module : G4hntools
# Package: Geant4.src.G4analysis.G4hntools
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4hntools
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
    g4hntools_defs.hh
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
)

# List any source specific properties here
