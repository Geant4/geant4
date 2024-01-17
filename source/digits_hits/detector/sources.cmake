# - G4detector module build definition

# Define the Geant4 Module.
geant4_add_module(G4detector
  PUBLIC_HEADERS
    G4CellScoreComposer.hh
    G4CellScoreValues.hh
    G4CollectionNameVector.hh
    G4HCtable.hh
    G4MultiFunctionalDetector.hh
    G4SDManager.hh
    G4SDStructure.hh
    G4SDmessenger.hh
    G4SensitiveVolumeList.hh
    G4TrackLogger.hh
    G4TScoreHistFiller.hh
    G4TScoreHistFiller.icc
    G4VPrimitivePlotter.hh
    G4VPrimitiveScorer.hh
    G4VReadOutGeometry.hh
    G4VSDFilter.hh
    G4VScoreHistFiller.hh
    G4VSensitiveDetector.hh
    G4MultiSensitiveDetector.hh
  SOURCES
    G4CellScoreComposer.cc
    G4HCtable.cc
    G4MultiFunctionalDetector.cc
    G4SDManager.cc
    G4SDStructure.cc
    G4SDmessenger.cc
    G4SensitiveVolumeList.cc
    G4TrackLogger.cc
    G4VPrimitiveScorer.cc
    G4VReadOutGeometry.cc
    G4VSDFilter.cc
    G4VScoreHistFiller.cc
    G4VSensitiveDetector.cc
    G4MultiSensitiveDetector.cc)

geant4_module_link_libraries(G4detector
  PUBLIC G4track G4hits G4geometrymng G4intercoms G4globman
  PRIVATE G4navigation)
