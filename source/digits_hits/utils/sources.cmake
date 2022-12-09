# - G4detutils module build definition

# Define the Geant4 Module.
geant4_add_module(G4detutils
  PUBLIC_HEADERS
    G4DefaultLinearColorMap.hh
    G4ScoreLogColorMap.hh
    G4VScoreNtupleWriter.hh
    G4TScoreNtupleWriter.hh
    G4TScoreNtupleWriter.icc
    G4TScoreNtupleWriterMessenger.hh
    G4TScoreNtupleWriterMessenger.icc
    G4ScoreQuantityMessenger.hh
    G4ScoringBox.hh
    G4ScoringCylinder.hh
    G4ScoringManager.hh
    G4ScoringMessenger.hh
    G4ScoringRealWorld.hh
    G4ScoringProbe.hh
    G4VScoreColorMap.hh
    G4VScoreWriter.hh
    G4VScoringMesh.hh
  SOURCES
    G4DefaultLinearColorMap.cc
    G4ScoreLogColorMap.cc
    G4VScoreNtupleWriter.cc
    G4ScoreQuantityMessenger.cc
    G4ScoringBox.cc
    G4ScoringCylinder.cc
    G4ScoringManager.cc
    G4ScoringMessenger.cc
    G4ScoringRealWorld.cc
    G4ScoringProbe.cc
    G4VScoreColorMap.cc
    G4VScoreWriter.cc
    G4VScoringMesh.cc)

geant4_module_link_libraries(G4detutils
  PUBLIC G4hits G4intercoms G4hepnumerics G4globman
  PRIVATE G4detector G4detscorer G4volumes G4geomdivision G4csg G4geometrymng G4graphics_reps G4materials)
