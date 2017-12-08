#------------------------------------------------------------------------------
# sources.cmake
# Module : G4emutils
# Package: Geant4.src.G4processes.G4electromagnetic.G4emutils
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 105801 2017-08-21 07:37:34Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4emutils
    HEADERS
        G4AngleDirect.hh
        G4AtomicShell.hh
        G4AtomicShellEnumerator.hh
        G4DummyModel.hh
        G4ElectronIonPair.hh
        G4EmBiasingManager.hh
        G4EmCalculator.hh
        G4EmConfigurator.hh
        G4EmCorrections.hh
        G4EmDataHandler.hh
        G4EmElementSelector.hh
        G4EmModelManager.hh
        G4EmMultiModel.hh
        G4EmParameters.hh
        G4EmParametersMessenger.hh
        G4EmProcessOptions.hh
        G4EmProcessSubType.hh
        G4EmSaturation.hh
        G4EmTableType.hh
        G4EnergyLossTables.hh
        G4LossTableBuilder.hh
        G4LossTableManager.hh
        G4MscStepLimitType.hh
        G4NuclearFormfactorType.hh
        G4VAtomDeexcitation.hh
        G4VEmAngularDistribution.hh
        G4VEmFluctuationModel.hh
        G4VEmModel.hh
        G4VEmProcess.hh
        G4VEnergyLossProcess.hh
        G4VMscModel.hh
        G4VMultipleScattering.hh
        G4VSubCutProducer.hh
        G4ionEffectiveCharge.hh
    SOURCES
        G4AngleDirect.cc
        G4DummyModel.cc
        G4ElectronIonPair.cc
        G4EmBiasingManager.cc
        G4EmCalculator.cc
        G4EmConfigurator.cc
        G4EmCorrections.cc
        G4EmDataHandler.cc
        G4EmElementSelector.cc 
        G4EmModelManager.cc
        G4EmMultiModel.cc
        G4EmParameters.cc
        G4EmParametersMessenger.cc
        G4EmProcessOptions.cc
        G4EmSaturation.cc
        G4EnergyLossTables.cc
        G4LossTableBuilder.cc
        G4LossTableManager.cc
        G4VAtomDeexcitation.cc
        G4VEmAngularDistribution.cc
        G4VEmFluctuationModel.cc
        G4VEmModel.cc
        G4VEmProcess.cc
        G4VEnergyLossProcess.cc
        G4VMscModel.cc
        G4VMultipleScattering.cc
        G4ionEffectiveCharge.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4cuts
        G4geometrymng
        G4globman
        G4intercoms
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4navigation
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

