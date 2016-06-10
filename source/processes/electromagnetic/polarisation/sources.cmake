#------------------------------------------------------------------------------
# sources.cmake
# Module : G4empolar
# Package: Geant4.src.G4processes.G4electromagnetic.G4empolar
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 66892 2013-01-17 10:57:59Z gunter $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/cuts/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/standard/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4empolar
    HEADERS
        G4PolarizationHelper.hh
        G4PolarizationManager.hh
        G4PolarizationMessenger.hh
        G4PolarizedAnnihilationCrossSection.hh
        G4PolarizedAnnihilationModel.hh
        G4PolarizedBhabhaCrossSection.hh
        G4PolarizedBremsstrahlungCrossSection.hh
        G4PolarizedCompton.hh
        G4PolarizedComptonCrossSection.hh
        G4PolarizedComptonModel.hh
        G4PolarizedGammaConversion.hh
        G4PolarizedGammaConversionModel.hh
        G4PolarizedMollerBhabhaModel.hh
        G4PolarizedMollerCrossSection.hh
        G4PolarizedPEEffectCrossSection.hh
        G4PolarizedPEEffectModel.hh
        G4PolarizedPairProductionCrossSection.hh
        G4PolarizedPhotoElectricEffect.hh
        G4StokesVector.hh
        G4VPolarizedCrossSection.hh
        G4ePolarizedBremsstrahlung.hh
        G4ePolarizedBremsstrahlungModel.hh
        G4ePolarizedIonisation.hh
        G4eplusPolarizedAnnihilation.hh
    SOURCES
        G4PolarizationHelper.cc
        G4PolarizationManager.cc
        G4PolarizationMessenger.cc
        G4PolarizedAnnihilationCrossSection.cc
        G4PolarizedAnnihilationModel.cc
        G4PolarizedBhabhaCrossSection.cc
        G4PolarizedBremsstrahlungCrossSection.cc
        G4PolarizedCompton.cc
        G4PolarizedComptonCrossSection.cc
        G4PolarizedComptonModel.cc
        G4PolarizedGammaConversion.cc
        G4PolarizedGammaConversionModel.cc
        G4PolarizedMollerBhabhaModel.cc
        G4PolarizedMollerCrossSection.cc
        G4PolarizedPEEffectCrossSection.cc
        G4PolarizedPEEffectModel.cc
        G4PolarizedPairProductionCrossSection.cc
        G4PolarizedPhotoElectricEffect.cc
        G4StokesVector.cc
        G4VPolarizedCrossSection.cc
        G4ePolarizedBremsstrahlung.cc
        G4ePolarizedBremsstrahlungModel.cc
        G4ePolarizedIonisation.cc
        G4eplusPolarizedAnnihilation.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4cuts
        G4emstandard
        G4emutils
        G4geometrymng
        G4globman
        G4hepnumerics
        G4intercoms
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4partman
        G4procman
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

