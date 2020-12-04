#------------------------------------------------------------------------------
# Module : G4detutils
# Package: Geant4.src.G4digits_hits.G4detutils
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4detutils
  HEADERS
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
    G4ScoreQuantityMessengerQCmd.cc
    G4ScoringBox.cc
    G4ScoringCylinder.cc
    G4ScoringManager.cc
    G4ScoringMessenger.cc
    G4ScoringRealWorld.cc
    G4ScoringProbe.cc
    G4VScoreColorMap.cc
    G4VScoreWriter.cc
    G4VScoringMesh.cc
  GRANULAR_DEPENDENCIES
    G4csg
    G4detector
    G4detscorer
    G4digits
    G4geomdivision
    G4geometrymng
    G4globman
    G4graphics_reps
    G4hits
    G4intercoms
    G4materials
    G4navigation
    G4partman
    G4track
    G4volumes
  GLOBAL_DEPENDENCIES
    G4geometry
    G4global
    G4graphics_reps
    G4intercoms
    G4materials
    G4particles
    G4track
)

# List any source specific properties here
