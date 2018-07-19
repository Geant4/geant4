#------------------------------------------------------------------------------
# sources.cmake
# Module : G4event
# Package: Geant4.src.G4event
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 95248 2016-02-02 09:37:04Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/detector/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/digits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/biasing/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPNumerics/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/barions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/ions/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/hadrons/mesons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/leptons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/tracking/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4event
    HEADERS
        G4AdjointPosOnPhysVolGenerator.hh
        G4AdjointPrimaryGenerator.hh
        G4AdjointStackingAction.hh
        G4ClassificationOfNewTrack.hh
        G4EvManMessenger.hh
        G4Event.hh
        G4EventManager.hh
        G4GeneralParticleSource.hh
        G4GeneralParticleSourceMessenger.hh
        G4HEPEvtInterface.hh
        G4HEPEvtParticle.hh
        G4ParticleGun.hh
        G4ParticleGunMessenger.hh
        G4PrimaryTransformer.hh
        G4RayShooter.hh
        G4SPSAngDistribution.hh
        G4SPSEneDistribution.hh
        G4SPSPosDistribution.hh
        G4SPSRandomGenerator.hh
        G4SingleParticleSource.hh
        G4SmartTrackStack.hh
        G4StackChecker.hh
        G4StackManager.hh
        G4StackedTrack.hh
        G4StackingMessenger.hh
        G4TrackStack.hh
        G4TrajectoryContainer.hh
        G4UserEventAction.hh
        G4MultiEventAction.hh
        G4UserStackingAction.hh
        G4VPrimaryGenerator.hh
        G4VUserEventInformation.hh
        eventgendefs.hh
        evmandefs.hh
        evtdefs.hh
        trajectoryControl.hh
	G4GeneralParticleSourceData.hh
    SOURCES
        G4AdjointPosOnPhysVolGenerator.cc
        G4AdjointPrimaryGenerator.cc
        G4AdjointStackingAction.cc
        G4EvManMessenger.cc
        G4Event.cc
        G4EventManager.cc
        G4GeneralParticleSource.cc
        G4GeneralParticleSourceMessenger.cc
        G4HEPEvtInterface.cc
        G4HEPEvtParticle.cc
        G4ParticleGun.cc
        G4ParticleGunMessenger.cc
        G4PrimaryTransformer.cc
        G4RayShooter.cc
        G4SPSAngDistribution.cc
        G4SPSEneDistribution.cc
        G4SPSPosDistribution.cc
        G4SPSRandomGenerator.cc
        G4SingleParticleSource.cc
        G4SmartTrackStack.cc
        G4StackChecker.cc
        G4StackManager.cc
        G4StackingMessenger.cc
        G4TrackStack.cc
        G4TrajectoryContainer.cc
        G4UserEventAction.cc
        G4MultiEventAction.cc
        G4UserStackingAction.cc
        G4VPrimaryGenerator.cc
	G4GeneralParticleSourceData.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4detector
        G4digits
        G4emutils
        G4geombias
        G4geometrymng
        G4globman
        G4graphics_reps
        G4hepnumerics
        G4hits
        G4intercoms
        G4ions
        G4leptons
        G4materials
        G4mesons
        G4navigation
        G4partman
        G4procman
        G4track
        G4tracking
        G4volumes
    GLOBAL_DEPENDENCIES
        G4digits_hits
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

