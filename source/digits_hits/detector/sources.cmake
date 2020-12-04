#------------------------------------------------------------------------------
# Module : G4detector
# Package: Geant4.src.G4digits_hits.G4detector
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4detector
  HEADERS
    G4CellScoreComposer.hh
    G4CellScoreValues.hh
    G4CollectionNameVector.hh
    G4HCtable.hh
    G4MultiFunctionalDetector.hh
    G4SDManager.hh
    G4SDStructure.hh
    G4SDmessenger.hh
    G4SensitiveVolumeList.hh
    G4SensitiveVolumeList.icc
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
    G4MultiSensitiveDetector.cc
  GRANULAR_DEPENDENCIES
    G4geombias
    G4geometrymng
    G4globman
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
    G4intercoms
    G4materials
    G4particles
    G4track
)

# List any source specific properties here
