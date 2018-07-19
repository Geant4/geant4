#------------------------------------------------------------------------------
# sources.cmake
# Module : G4emadjoint
# Package: Geant4.src.G4processes.G4electromagnetic.G4emadjoint
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 100666 2016-10-31 10:27:00Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/adjoint/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
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
GEANT4_DEFINE_MODULE(NAME G4emadjoint
    HEADERS
        G4AdjointAlongStepWeightCorrection.hh
        G4AdjointBremsstrahlungModel.hh
        G4AdjointCSManager.hh
        G4AdjointCSMatrix.hh
        G4AdjointComptonModel.hh
        G4AdjointInterpolator.hh
        G4AdjointIonIonisationModel.hh
        G4AdjointPhotoElectricModel.hh
        G4AdjointProcessEquivalentToDirectProcess.hh
        G4AdjointeIonisationModel.hh
        G4AdjointhIonisationModel.hh
        G4AdjointhMultipleScattering.hh
        G4ContinuousGainOfEnergy.hh
        G4InversePEEffect.hh
        G4IonInverseIonisation.hh
        G4VAdjointReverseReaction.hh
        G4AdjointForcedInteractionForGamma.hh
        G4VEmAdjointModel.hh
        G4eInverseBremsstrahlung.hh
        G4eInverseCompton.hh
        G4eInverseIonisation.hh
        G4hInverseIonisation.hh
        G4UrbanAdjointMscModel.hh
        G4eAdjointMultipleScattering.hh
    SOURCES
        G4AdjointAlongStepWeightCorrection.cc
        G4AdjointBremsstrahlungModel.cc
        G4AdjointCSManager.cc
        G4AdjointCSMatrix.cc
        G4AdjointComptonModel.cc
        G4AdjointInterpolator.cc
        G4AdjointIonIonisationModel.cc
        G4AdjointPhotoElectricModel.cc
        G4AdjointProcessEquivalentToDirectProcess.cc
        G4AdjointeIonisationModel.cc
        G4AdjointhIonisationModel.cc
        G4AdjointhMultipleScattering.cc
        G4ContinuousGainOfEnergy.cc
        G4InversePEEffect.cc
        G4IonInverseIonisation.cc
        G4VAdjointReverseReaction.cc
        G4AdjointForcedInteractionForGamma.cc
        G4VEmAdjointModel.cc
        G4eInverseBremsstrahlung.cc
        G4eInverseCompton.cc
        G4eInverseIonisation.cc
        G4hInverseIonisation.cc
        G4UrbanAdjointMscModel.cc
        G4eAdjointMultipleScattering.cc
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
        G4navigation
        G4partadj
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

