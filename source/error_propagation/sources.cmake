#------------------------------------------------------------------------------
# sources.cmake
# Module : G4error_propagation
# Package: Geant4.src.G4error_propagation
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
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/digits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/event/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/magneticfield/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/solids/CSG/include)
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
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/muons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/standard/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/electromagnetic/utils/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/transportation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/run/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/tracking/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4error_propagation
    HEADERS
        G4ErrorFreeTrajParam.hh
        G4ErrorFreeTrajState.hh
        G4ErrorGeomVolumeTarget.hh
        G4ErrorMagFieldLimitProcess.hh
        G4ErrorMatrix.hh
        G4ErrorMatrix.icc
        G4ErrorMessenger.hh
        G4ErrorPhysicsList.hh
        G4ErrorPropagator.hh
        G4ErrorPropagatorManager.hh
        G4ErrorRunManagerHelper.hh
        G4ErrorStepLengthLimitProcess.hh
        G4ErrorSurfaceTrajParam.hh
        G4ErrorSurfaceTrajState.hh
        G4ErrorSymMatrix.hh
        G4ErrorSymMatrix.icc
        G4ErrorTrackLengthTarget.hh
        G4ErrorTrajErr.hh
        G4ErrorTrajState.hh
        G4VErrorLimitProcess.hh
    SOURCES
        G4ErrorFreeTrajParam.cc
        G4ErrorFreeTrajState.cc
        G4ErrorGeomVolumeTarget.cc
        G4ErrorMagFieldLimitProcess.cc
        G4ErrorMatrix.cc
        G4ErrorMessenger.cc
        G4ErrorPhysicsList.cc
        G4ErrorPropagator.cc
        G4ErrorPropagatorManager.cc
        G4ErrorRunManagerHelper.cc
        G4ErrorStepLengthLimitProcess.cc
        G4ErrorSurfaceTrajParam.cc
        G4ErrorSurfaceTrajState.cc
        G4ErrorSymMatrix.cc
        G4ErrorTrackLengthTarget.cc
        G4ErrorTrajState.cc
        G4VErrorLimitProcess.cc
    GRANULAR_DEPENDENCIES
        G4baryons
        G4bosons
        G4csg
        G4cuts
        G4digits
        G4emstandard
        G4emutils
        G4event
        G4geometrymng
        G4globman
        G4hits
        G4intercoms
        G4ions
        G4leptons
        G4magneticfield
        G4materials
        G4mesons
        G4muons
        G4navigation
        G4partman
        G4procman
        G4run
        G4track
        G4tracking
        G4transportation
        G4volumes
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4event
        G4geometry
        G4global
        G4intercoms
        G4materials
        G4particles
        G4processes
        G4run
        G4track
        G4tracking
    LINK_LIBRARIES
)

# List any source specific properties here

