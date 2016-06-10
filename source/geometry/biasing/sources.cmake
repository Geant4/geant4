#------------------------------------------------------------------------------
# sources.cmake
# Module : G4geombias
# Package: Geant4.src.G4geometry.G4geombias
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 66356 2012-12-18 09:02:32Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4geombias
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
    LINK_LIBRARIES
)

# List any source specific properties here

