#------------------------------------------------------------------------------
# sources.cmake
# Module : G4materials
# Package: Geant4.src.G4materials
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 106115 2017-09-13 10:16:51Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4materials
    HEADERS
			G4AtomicBond.hh
			G4AtomicFormFactor.hh
			G4AtomicShells.hh
			G4AtomicShells_EADL.hh 
			G4AtomicShells_EADL.hh 
			G4AtomicShells_EADL.hh 
			G4CrystalAtomBase.hh
			G4CrystalBravaisLattices.h
			G4CrystalLatticeSystems.h
			G4CrystalExtension.hh
			G4CrystalUnitCell.hh
			G4DensityEffectData.hh
			G4ElementData.hh
			G4Element.hh
			G4ElementTable.hh
			G4ElementVector.hh
			G4ExtDEDXTable.hh
			G4ExtendedMaterial.hh
			G4IonisParamElm.hh
			G4IonisParamMat.hh
			G4IonStoppingData.hh
			G4Isotope.hh
			G4IsotopeVector.hh
			G4LatticeLogical.hh
			G4LatticePhysical.hh
			G4Material.hh
                        G4MaterialPropertiesIndex.hh
			G4MaterialPropertiesTable.hh
			G4MaterialPropertiesTable.icc
			G4MaterialPropertyVector.hh
			G4MaterialTable.hh
			G4NistElementBuilder.hh
			G4NistManager.hh
			G4NistMaterialBuilder.hh
			G4NistMessenger.hh
			G4OpticalSurface.hh
			G4SandiaTable.hh
			G4StaticSandiaData.hh
			G4SurfaceProperty.hh
			G4UCNMaterialPropertiesTable.hh
			G4UCNMicroRoughnessHelper.hh
			G4VIonDEDXTable.hh
			G4VMaterialExtension.hh
    SOURCES
			G4AtomicBond.cc
			G4AtomicFormFactor.cc
			G4AtomicShells.cc
			G4AtomicShells_EADL.cc 
			G4AtomicShells_EADL.cc 
			G4AtomicShells_EADL.cc 
			G4CrystalExtension.cc
			G4CrystalUnitCell.cc
			G4DensityEffectData.cc
			G4Element.cc
			G4ElementData.cc
			G4ExtDEDXTable.cc
			G4ExtendedMaterial.cc
			G4IonisParamElm.cc
			G4IonisParamMat.cc
			G4IonStoppingData.cc
			G4Isotope.cc
			G4LatticeLogical.cc
			G4LatticePhysical.cc
			G4Material.cc
			G4MaterialPropertiesTable.cc
			G4NistElementBuilder.cc
			G4NistManager.cc
			G4NistMaterialBuilder.cc
			G4NistMessenger.cc
			G4OpticalSurface.cc
			G4SandiaTable.cc
			G4SurfaceProperty.cc
			G4UCNMaterialPropertiesTable.cc
			G4UCNMicroRoughnessHelper.cc
			G4VIonDEDXTable.cc
    GRANULAR_DEPENDENCIES
        G4globman
        G4intercoms
    GLOBAL_DEPENDENCIES
        G4global
        G4intercoms
    LINK_LIBRARIES
)

# List any source specific properties here
