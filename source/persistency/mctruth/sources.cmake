#------------------------------------------------------------------------------
# sources.cmake
# Module : G4mctruth
# Package: Geant4.src.G4persistency.G4mctruth
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
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/run/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/tracking/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4mctruth
    HEADERS
        G4DCIOcatalog.hh
        G4DCIOentryT.hh
        G4FileUtilities.hh
        G4HCIOcatalog.hh
        G4HCIOentryT.hh
        G4MCTEvent.hh
        G4MCTGenEvent.hh
        G4MCTGenParticle.hh
        G4MCTSimEvent.hh
        G4MCTSimParticle.hh
        G4MCTSimVertex.hh
        G4PersistencyCenter.hh
        G4PersistencyCenterMessenger.hh
        G4PersistencyManager.hh
        G4PersistencyManagerT.hh
        G4Pevent.hh
        G4VDCIOentry.hh
        G4VHCIOentry.hh
        G4VMCTruthIO.hh
        G4VPDigitIO.hh
        G4VPDigitsCollectionIO.hh
        G4VPEventIO.hh
        G4VPHitIO.hh
        G4VPHitsCollectionIO.hh
        G4VTransactionManager.hh
    SOURCES
        G4DCIOcatalog.cc
        G4FileUtilities.cc
        G4HCIOcatalog.cc
        G4MCTEvent.cc
        G4MCTGenEvent.cc
        G4MCTSimEvent.cc
        G4MCTSimParticle.cc
        G4MCTSimVertex.cc
        G4PersistencyCenter.cc
        G4PersistencyCenterMessenger.cc
        G4PersistencyManager.cc
        G4Pevent.cc
        G4VDCIOentry.cc
        G4VHCIOentry.cc
        G4VMCTruthIO.cc
        G4VPDigitIO.cc
        G4VPDigitsCollectionIO.cc
        G4VPEventIO.cc
        G4VPHitIO.cc
        G4VPHitsCollectionIO.cc
    GRANULAR_DEPENDENCIES
        G4digits
        G4event
        G4geometrymng
        G4globman
        G4hits
        G4intercoms
        G4partman
        G4run
        G4track
        G4tracking
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4event
        G4geometry
        G4global
        G4intercoms
        G4particles
        G4run
        G4track
        G4tracking
    LINK_LIBRARIES
)

# List any source specific properties here

