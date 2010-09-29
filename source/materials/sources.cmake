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
# $Id: sources.cmake,v 1.1 2010-09-29 18:46:35 bmorgan Exp $
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
        G4AtomicShells.hh
        G4DensityEffectData.hh
        G4Element.hh
        G4ElementTable.hh
        G4ElementVector.hh
        G4ExtDEDXTable.hh
        G4IonStoppingData.hh
        G4IonisParamElm.hh
        G4IonisParamMat.hh
        G4IronStoppingICRU73.hh
        G4Isotope.hh
        G4IsotopeVector.hh
        G4MPVEntry.hh
        G4Material.hh
        G4MaterialPropertiesTable.hh
        G4MaterialPropertiesTable.icc
        G4MaterialPropertyVector.hh
        G4MaterialPropertyVector.icc
        G4MaterialStoppingICRU73.hh
        G4MaterialTable.hh
        G4NistElementBuilder.hh
        G4NistManager.hh
        G4NistMaterialBuilder.hh
        G4NistMessenger.hh
        G4OpticalSurface.hh
        G4SandiaTable.hh
        G4SimpleMaterialStoppingICRU73.hh
        G4StaticSandiaData.hh
        G4SurfaceProperty.hh
        G4VIonDEDXTable.hh
    SOURCES
        G4AtomicShells.cc
        G4DensityEffectData.cc
        G4Element.cc
        G4ExtDEDXTable.cc
        G4IonStoppingData.cc
        G4IonisParamElm.cc
        G4IonisParamMat.cc
        G4IronStoppingICRU73.cc
        G4Isotope.cc
        G4MPVEntry.cc
        G4Material.cc
        G4MaterialPropertiesTable.cc
        G4MaterialPropertyVector.cc
        G4MaterialStoppingICRU73.cc
        G4NistElementBuilder.cc
        G4NistManager.cc
        G4NistMaterialBuilder.cc
        G4NistMessenger.cc
        G4OpticalSurface.cc
        G4SandiaTable.cc
        G4SimpleMaterialStoppingICRU73.cc
        G4SurfaceProperty.cc
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

