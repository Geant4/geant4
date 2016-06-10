#------------------------------------------------------------------------------
# sources.cmake
# Module : G4geomtext
# Package: Geant4.src.G4persistency.G4geomtext
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 66363 2012-12-18 09:12:54Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/divisions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/Boolean/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/specific/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4geomtext
    HEADERS
        G4tgbDetectorBuilder.hh
        G4tgbDetectorConstruction.hh
        G4tgbElement.hh
        G4tgbGeometryDumper.hh
        G4tgbIsotope.hh
        G4tgbMaterial.hh
        G4tgbMaterialMgr.hh
        G4tgbMaterialMixture.hh
        G4tgbMaterialMixtureByNoAtoms.hh
        G4tgbMaterialMixtureByVolume.hh
        G4tgbMaterialMixtureByWeight.hh
        G4tgbMaterialSimple.hh
        G4tgbPlaceParamCircle.hh
        G4tgbPlaceParamLinear.hh
        G4tgbPlaceParamSquare.hh
        G4tgbPlaceParameterisation.hh
        G4tgbRotationMatrix.hh
        G4tgbRotationMatrixMgr.hh
        G4tgbVolume.hh
        G4tgbVolumeMgr.hh
        G4tgrElement.hh
        G4tgrElementFromIsotopes.hh
        G4tgrElementSimple.hh
        G4tgrEvaluator.hh
        G4tgrFileIn.hh
        G4tgrFileReader.hh
        G4tgrIsotope.hh
        G4tgrLineProcessor.hh
        G4tgrMaterial.hh
        G4tgrMaterialFactory.hh
        G4tgrMaterialMixture.hh
        G4tgrMaterialSimple.hh
        G4tgrMessenger.hh
        G4tgrParameterMgr.hh
        G4tgrPlace.hh
        G4tgrPlaceDivRep.hh
        G4tgrPlaceParameterisation.hh
        G4tgrPlaceSimple.hh
        G4tgrRotationMatrix.hh
        G4tgrRotationMatrixFactory.hh
        G4tgrSolid.hh
        G4tgrSolidBoolean.hh
        G4tgrUtils.hh
        G4tgrVolume.hh
        G4tgrVolumeAssembly.hh
        G4tgrVolumeDivision.hh
        G4tgrVolumeMgr.hh
    SOURCES
        G4tgbDetectorBuilder.cc
        G4tgbDetectorConstruction.cc
        G4tgbElement.cc
        G4tgbGeometryDumper.cc
        G4tgbIsotope.cc
        G4tgbMaterial.cc
        G4tgbMaterialMgr.cc
        G4tgbMaterialMixture.cc
        G4tgbMaterialMixtureByNoAtoms.cc
        G4tgbMaterialMixtureByVolume.cc
        G4tgbMaterialMixtureByWeight.cc
        G4tgbMaterialSimple.cc
        G4tgbPlaceParamCircle.cc
        G4tgbPlaceParamLinear.cc
        G4tgbPlaceParamSquare.cc
        G4tgbPlaceParameterisation.cc
        G4tgbRotationMatrix.cc
        G4tgbRotationMatrixMgr.cc
        G4tgbVolume.cc
        G4tgbVolumeMgr.cc
        G4tgrElement.cc
        G4tgrElementFromIsotopes.cc
        G4tgrElementSimple.cc
        G4tgrEvaluator.cc
        G4tgrFileIn.cc
        G4tgrFileReader.cc
        G4tgrIsotope.cc
        G4tgrLineProcessor.cc
        G4tgrMaterial.cc
        G4tgrMaterialFactory.cc
        G4tgrMaterialMixture.cc
        G4tgrMaterialSimple.cc
        G4tgrMessenger.cc
        G4tgrParameterMgr.cc
        G4tgrPlace.cc
        G4tgrPlaceDivRep.cc
        G4tgrPlaceParameterisation.cc
        G4tgrPlaceSimple.cc
        G4tgrRotationMatrix.cc
        G4tgrRotationMatrixFactory.cc
        G4tgrSolid.cc
        G4tgrSolidBoolean.cc
        G4tgrUtils.cc
        G4tgrVolume.cc
        G4tgrVolumeAssembly.cc
        G4tgrVolumeDivision.cc
        G4tgrVolumeMgr.cc
    GRANULAR_DEPENDENCIES
        G4csg
        G4geomBoolean
        G4geomdivision
        G4geometrymng
        G4globman
        G4graphics_reps
        G4hepnumerics
        G4intercoms
        G4materials
        G4partman
        G4specsolids
        G4volumes
    GLOBAL_DEPENDENCIES
        G4geometry
        G4global
        G4graphics_reps
        G4intercoms
        G4materials
        G4particles
    LINK_LIBRARIES
)

# List any source specific properties here

