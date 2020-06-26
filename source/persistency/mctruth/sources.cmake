#------------------------------------------------------------------------------
# Module : G4mctruth
# Package: Geant4.src.G4persistency.G4mctruth
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4mctruth
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

