#------------------------------------------------------------------------------
# sources.cmake
# Module : G4detscorer
# Package: Geant4.src.G4digits_hits.G4detscorer
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 88188 2015-02-02 17:22:17Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${USOLIDS_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/detector/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/digits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
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
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4detscorer
    HEADERS
        G4PSCellCharge.hh
        G4PSCellCharge3D.hh
        G4PSCellFlux.hh
        G4PSCellFlux3D.hh
        G4PSCellFluxForCylinder3D.hh
        G4PSCylinderSurfaceCurrent.hh
        G4PSCylinderSurfaceCurrent3D.hh
        G4PSCylinderSurfaceFlux.hh
        G4PSCylinderSurfaceFlux3D.hh
        G4PSDirectionFlag.hh
        G4PSDoseDeposit.hh
        G4PSDoseDeposit3D.hh
        G4PSDoseDepositForCylinder3D.hh
        G4PSEnergyDeposit.hh
        G4PSEnergyDeposit3D.hh
        G4PSFlatSurfaceCurrent.hh
        G4PSFlatSurfaceCurrent3D.hh
        G4PSFlatSurfaceFlux.hh
        G4PSFlatSurfaceFlux3D.hh
        G4PSMinKinEAtGeneration.hh
        G4PSMinKinEAtGeneration3D.hh
        G4PSNofCollision.hh
        G4PSNofCollision3D.hh
        G4PSNofSecondary.hh
        G4PSNofSecondary3D.hh
        G4PSNofStep.hh
        G4PSNofStep3D.hh
        G4PSPassageCellCurrent.hh
        G4PSPassageCellCurrent3D.hh
        G4PSPassageCellFlux.hh
        G4PSPassageCellFlux3D.hh
        G4PSPassageCellFluxForCylinder3D.hh
        G4PSPassageTrackLength.hh
        G4PSPassageTrackLength3D.hh
        G4PSPopulation.hh
        G4PSPopulation3D.hh
        G4PSSphereSurfaceCurrent.hh
        G4PSSphereSurfaceCurrent3D.hh
        G4PSSphereSurfaceFlux.hh
        G4PSSphereSurfaceFlux3D.hh
        G4PSStepChecker.hh
        G4PSStepChecker3D.hh
        G4PSTermination.hh
        G4PSTermination3D.hh
        G4PSTrackCounter.hh
        G4PSTrackCounter3D.hh
        G4PSTrackLength.hh
        G4PSTrackLength3D.hh
        G4SDChargedFilter.hh
        G4SDKineticEnergyFilter.hh
        G4SDNeutralFilter.hh
        G4SDParticleFilter.hh
        G4SDParticleWithEnergyFilter.hh
    SOURCES
        G4PSCellCharge.cc
        G4PSCellCharge3D.cc
        G4PSCellFlux.cc
        G4PSCellFlux3D.cc
        G4PSCellFluxForCylinder3D.cc
        G4PSCylinderSurfaceCurrent.cc
        G4PSCylinderSurfaceCurrent3D.cc
        G4PSCylinderSurfaceFlux.cc
        G4PSCylinderSurfaceFlux3D.cc
        G4PSDoseDeposit.cc
        G4PSDoseDeposit3D.cc
        G4PSDoseDepositForCylinder3D.cc
        G4PSEnergyDeposit.cc
        G4PSEnergyDeposit3D.cc
        G4PSFlatSurfaceCurrent.cc
        G4PSFlatSurfaceCurrent3D.cc
        G4PSFlatSurfaceFlux.cc
        G4PSFlatSurfaceFlux3D.cc
        G4PSMinKinEAtGeneration.cc
        G4PSMinKinEAtGeneration3D.cc
        G4PSNofCollision.cc
        G4PSNofCollision3D.cc
        G4PSNofSecondary.cc
        G4PSNofSecondary3D.cc
        G4PSNofStep.cc
        G4PSNofStep3D.cc
        G4PSPassageCellCurrent.cc
        G4PSPassageCellCurrent3D.cc
        G4PSPassageCellFlux.cc
        G4PSPassageCellFlux3D.cc
        G4PSPassageCellFluxForCylinder3D.cc
        G4PSPassageTrackLength.cc
        G4PSPassageTrackLength3D.cc
        G4PSPopulation.cc
        G4PSPopulation3D.cc
        G4PSSphereSurfaceCurrent.cc
        G4PSSphereSurfaceCurrent3D.cc
        G4PSSphereSurfaceFlux.cc
        G4PSSphereSurfaceFlux3D.cc
        G4PSStepChecker.cc
        G4PSStepChecker3D.cc
        G4PSTermination.cc
        G4PSTermination3D.cc
        G4PSTrackCounter.cc
        G4PSTrackCounter3D.cc
        G4PSTrackLength.cc
        G4PSTrackLength3D.cc
        G4SDChargedFilter.cc
        G4SDKineticEnergyFilter.cc
        G4SDNeutralFilter.cc
        G4SDParticleFilter.cc
        G4SDParticleWithEnergyFilter.cc
    GRANULAR_DEPENDENCIES
        G4csg
        G4detector
        G4digits
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
    LINK_LIBRARIES
)

# List any source specific properties here

