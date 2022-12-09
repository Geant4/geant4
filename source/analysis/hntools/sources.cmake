# - G4hntools module build definition

# Define the Geant4 Module.
geant4_add_module(G4hntools
  PUBLIC_HEADERS
    G4MPIToolsManager.hh
    G4MPIToolsManager.icc
    G4PlotManager.hh
    G4PlotMessenger.hh
    G4PlotParameters.hh
    G4TH1ToolsManager.hh
    G4TH2ToolsManager.hh
    G4TH3ToolsManager.hh
    G4TP1ToolsManager.hh
    G4TP2ToolsManager.hh
    G4ToolsAnalysisManager.hh
    G4ToolsAnalysisManager.icc
    G4ToolsAnalysisReader.hh
    G4ToolsAnalysisReader.icc
    g4hntools_defs.hh
  SOURCES
    G4PlotManager.cc
    G4PlotMessenger.cc
    G4PlotParameters.cc
    G4TH1ToolsManager.cc
    G4TH2ToolsManager.cc
    G4TH3ToolsManager.cc
    G4TP1ToolsManager.cc
    G4TP2ToolsManager.cc
    G4ToolsAnalysisManager.cc
    G4ToolsAnalysisReader.cc)

geant4_module_link_libraries(G4hntools PUBLIC G4analysismng G4globman G4intercoms G4tools)
