# - G4scoring module build definition

# Define the Geant4 Module.
geant4_add_module(G4scoring
  PUBLIC_HEADERS
    G4EnergySplitter.hh
    G4EnergySplitter.icc
    G4ParallelWorldProcess.hh
    G4ParallelWorldProcessStore.hh
    G4ParallelWorldScoringProcess.hh
    G4ScoreSplittingProcess.hh
  SOURCES
    G4EnergySplitter.cc
    G4ParallelWorldProcess.cc
    G4ParallelWorldProcessStore.cc
    G4ParallelWorldScoringProcess.cc
    G4ScoreSplittingProcess.cc)

geant4_module_link_libraries(G4scoring
  PUBLIC
    G4geometrymng
    G4globman
    G4intercoms
    G4magneticfield
    G4navigation
    G4procman
  PRIVATE
    G4cuts
    G4detector
    G4emutils
    G4materials
    G4muons
    G4partman
    G4track
    G4transportation
    G4volumes)
