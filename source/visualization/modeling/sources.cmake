#------------------------------------------------------------------------------
# sources.cmake
# Module : G4modeling
# Package: Geant4.src.G4visualization.G4modeling
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 103926 2017-05-03 13:43:27Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${USOLIDS_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/detector/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/digits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/event/include)
include_directories(${CMAKE_SOURCE_DIR}/source/run/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/magneticfield/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/Boolean/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/tracking/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4modeling
    HEADERS
        G4ArrowModel.hh
        G4AttFilterUtils.hh
        G4AttValueFilterT.hh
        G4AttributeFilterT.hh
        G4AxesModel.hh
        G4BoundingSphereScene.hh
        G4CallbackModel.hh
        G4DigiFilterFactories.hh
        G4DigiModel.hh
        G4GPSModel.hh
        G4HitFilterFactories.hh
        G4HitsModel.hh
        G4LogicalVolumeModel.hh
        G4MagneticFieldModel.hh
        G4ModelApplyCommandsT.hh
        G4ModelColourMap.hh
        G4ModelCommandUtils.hh
        G4ModelCommandsT.hh
        G4ModelCompoundCommandsT.hh
        G4ModelingParameters.hh
        G4ModelingParameters.icc
        G4NullModel.hh
        G4PSHitsModel.hh
        G4PhysicalVolumeMassScene.hh
        G4PhysicalVolumeModel.hh
        G4PhysicalVolumeSearchScene.hh
        G4PseudoScene.hh
        G4ScaleModel.hh
        G4TextModel.hh
        G4TouchableDumpScene.hh
        G4TrajectoriesModel.hh
        G4TrajectoryChargeFilter.hh
        G4TrajectoryDrawByAttribute.hh
        G4TrajectoryDrawByCharge.hh
        G4TrajectoryDrawByOriginVolume.hh
        G4TrajectoryDrawByParticleID.hh
        G4TrajectoryDrawByEncounteredVolume.hh
        G4TrajectoryDrawerUtils.hh
        G4TrajectoryFilterFactories.hh
        G4TrajectoryGenericDrawer.hh
        G4TrajectoryModelFactories.hh
        G4TrajectoryOriginVolumeFilter.hh
        G4TrajectoryParticleFilter.hh
        G4TrajectoryEncounteredVolumeFilter.hh
        G4VAttValueFilter.hh
        G4VModel.hh
        G4VModel.icc
        G4VModelCommand.hh
        G4VModelFactory.hh
        G4VTrajectoryModel.hh
        G4VisTrajContext.hh
        G4VisTrajContext.icc
    SOURCES
        G4ArrowModel.cc
        G4AttFilterUtils.cc
        G4AxesModel.cc
        G4BoundingSphereScene.cc
        G4DigiFilterFactories.cc
        G4DigiModel.cc
        G4GPSModel.cc
        G4HitFilterFactories.cc
        G4HitsModel.cc
        G4LogicalVolumeModel.cc
        G4MagneticFieldModel.cc
        G4ModelingParameters.cc
        G4NullModel.cc
        G4PSHitsModel.cc
        G4PhysicalVolumeMassScene.cc
        G4PhysicalVolumeModel.cc
        G4PhysicalVolumeSearchScene.cc
        G4ScaleModel.cc
        G4TextModel.cc
        G4TouchableDumpScene.cc
        G4TrajectoriesModel.cc
        G4TrajectoryChargeFilter.cc
        G4TrajectoryDrawByAttribute.cc
        G4TrajectoryDrawByCharge.cc
        G4TrajectoryDrawByOriginVolume.cc
        G4TrajectoryDrawByParticleID.cc
        G4TrajectoryDrawByEncounteredVolume.cc
        G4TrajectoryDrawerUtils.cc
        G4TrajectoryFilterFactories.cc
        G4TrajectoryGenericDrawer.cc
        G4TrajectoryModelFactories.cc
        G4TrajectoryOriginVolumeFilter.cc
        G4TrajectoryParticleFilter.cc
        G4TrajectoryEncounteredVolumeFilter.cc
        G4VModel.cc
        G4VTrajectoryModel.cc
        G4VisTrajContext.cc
    GRANULAR_DEPENDENCIES
        G4csg
        G4detector
        G4detutils
        G4digits
        G4run
        G4event
        G4geomBoolean
        G4geometrymng
        G4globman
        G4graphics_reps
        G4hepnumerics
        G4hits
        G4intercoms
        G4materials
        G4navigation
        G4partman
        G4procman
        G4specsolids
        G4track
        G4tracking
        G4volumes
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4run
        G4event
        G4geometry
        G4global
        G4graphics_reps
        G4intercoms
        G4materials
        G4particles
        G4processes
        G4track
        G4tracking
    LINK_LIBRARIES
)

# List any source specific properties here

