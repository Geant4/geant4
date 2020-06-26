#------------------------------------------------------------------------------
# Module : G4geombias
# Package: Geant4.src.G4geometry.G4geombias
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4geombias
  HEADERS
    G4GeometryCell.hh
    G4GeometryCellComp.hh
    G4GeometryCellImportance.hh
    G4GeometryCellStep.hh
    G4GeometryCellStepStream.hh
    G4GeometryCellWeight.hh
    G4IStore.hh
    G4ImportanceAlgorithm.hh
    G4Nsplit_Weight.hh
    G4VGCellFinder.hh
    G4VIStore.hh
    G4VImportanceAlgorithm.hh
    G4VImportanceSplitExaminer.hh
    G4VWeightWindowAlgorithm.hh
    G4VWeightWindowStore.hh
    G4WeightWindowAlgorithm.hh
    G4WeightWindowStore.hh
  SOURCES
    G4GeometryCell.cc
    G4GeometryCellComp.cc
    G4GeometryCellImportance.cc
    G4GeometryCellStep.cc
    G4GeometryCellStepStream.cc
    G4IStore.cc
    G4ImportanceAlgorithm.cc
    G4Nsplit_Weight.cc
    G4VGCellFinder.cc
    G4VIStore.cc
    G4VImportanceAlgorithm.cc
    G4VImportanceSplitExaminer.cc
    G4VWeightWindowAlgorithm.cc
    G4VWeightWindowStore.cc
    G4WeightWindowAlgorithm.cc
    G4WeightWindowStore.cc
  GRANULAR_DEPENDENCIES
    G4csg
    G4geometrymng
    G4globman
    G4graphics_reps
    G4intercoms
    G4materials
    G4navigation
    G4volumes
  GLOBAL_DEPENDENCIES
    G4global
    G4graphics_reps
    G4intercoms
    G4materials
)

# List any source specific properties here
