# - G4hntools module build definition

# Define the Geant4 Module.
geant4_add_module(G4hntools
  PUBLIC_HEADERS
    G4BaseHistoUtilities.hh
    G4MPIToolsManager.hh
    G4MPIToolsManager.icc
    G4H1ToolsManager.hh
    G4H2ToolsManager.hh
    G4H3ToolsManager.hh
    G4P1ToolsManager.hh
    G4P2ToolsManager.hh
    G4ToolsAnalysisManager.hh
    G4ToolsAnalysisManager.icc
    G4ToolsAnalysisMessenger.hh
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
    G4ToolsAnalysisMessenger.cc
    G4ToolsAnalysisReader.cc)

geant4_module_link_libraries(G4hntools PUBLIC G4analysismng G4globman G4tools)
