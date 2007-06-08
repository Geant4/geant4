#include "DetectorSliceEventAction.hh"

#include "G4Event.hh"
#include "G4ios.hh"

#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"

#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "DetectorSliceSensitiveEmCalo.hh"
#include "DetectorSliceSensitiveHadCalo.hh"
#include "DetectorSliceCalorimeterHit.hh"

#include "DetectorSliceSteppingAction.hh"
#include "G4RunManager.hh"
#include "G4Timer.hh"

#include "DetectorSliceAnalysis.hh"


DetectorSliceEventAction::DetectorSliceEventAction() : 
  theSteppingAction( 0 ) {
  instanciateSteppingAction();
  eventTimer = new G4Timer;
}


DetectorSliceEventAction::~DetectorSliceEventAction() {
  delete eventTimer;
}


void DetectorSliceEventAction::BeginOfEventAction( const G4Event* ) {
  //G4cout << "\n---> Begin of event: " << evt->GetEventID() << G4endl;
  eventTimer->Start();
}


void DetectorSliceEventAction::EndOfEventAction( const G4Event* evt ) {

  // Get the IDs of the collections of calorimeter hits. 
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  G4int calorimeterCollID = 0;

  // Extract the hits.
  G4HCofThisEvent* collectionEventHits = evt->GetHCofThisEvent();
  DetectorSliceCalorimeterHitsCollection* collectionCaloHits = 0;

  // Hits of the EM calorimeter.
  calorimeterCollID = SDman->GetCollectionID( colNam="emCalCollection" );
  G4int numEmCalorimeterHits = 0;
  G4double energyDepositedInActiveEmCalorimeterLayers = 0.0;
  if ( collectionEventHits ) {
    collectionCaloHits = dynamic_cast< DetectorSliceCalorimeterHitsCollection* >
      ( collectionEventHits->GetHC( calorimeterCollID ) );
  }
  if ( collectionCaloHits ) {
    numEmCalorimeterHits = collectionCaloHits->entries();
    for ( G4int i=0; i < numEmCalorimeterHits; i++ ) {
      energyDepositedInActiveEmCalorimeterLayers += (*collectionCaloHits)[i]->GetEdep(); 
      //G4cout << " DetectorSliceEventAction::EndOfEventAction : " << G4endl
      //       << " \t EM iHit=" << i 
      //       << "  energy=" << (*collectionCaloHits)[i]->GetEdep()
      //       << G4endl;           //***DEBUG***
    }
  }

  // Hits of the HAD calorimeter.
  calorimeterCollID = SDman->GetCollectionID( colNam="hadCalCollection" );
  G4int numHadCalorimeterHits = 0;
  G4double energyDepositedInActiveHadCalorimeterLayers = 0.0;
  if ( collectionEventHits ) {
    collectionCaloHits = dynamic_cast< DetectorSliceCalorimeterHitsCollection* >
      ( collectionEventHits->GetHC( calorimeterCollID ) );
  }
  if ( collectionCaloHits ) {
    numHadCalorimeterHits = collectionCaloHits->entries();
    for ( G4int i=0; i < numHadCalorimeterHits; i++ ) {
      energyDepositedInActiveHadCalorimeterLayers += (*collectionCaloHits)[i]->GetEdep();
      //G4cout << " DetectorSliceEventAction::EndOfEventAction : " << G4endl
      //       << " \t HAD iHit=" << i 
      //       << "  energy=" << (*collectionCaloHits)[i]->GetEdep()
      //       << G4endl;           //***DEBUG***
    }
  }

  //G4cout << " DetectorSliceEventAction::EndOfEventAction : " //***DEBUG***
  //       << " \t numEmCalorimeterHits  = " << numEmCalorimeterHits << G4endl
  //       << " \t numHadCalorimeterHits = " << numHadCalorimeterHits << G4endl; 

  // Get the information about the total release of energy in all
  // subdetectors (Tracker, EM Cal, HAD Cal, Muon detector),
  // the eventual exiting radius of a primary muon,
  // and the Id and kinetic energy of the incident (primary) particle.
  // All these information come from the Stepping Action.
  G4double totalEdepAllParticlesInTracker = 
    theSteppingAction->getTotalEdepAllParticlesInTracker();
  G4double totalEdepAllParticlesInEmCal = 
    theSteppingAction->getTotalEdepAllParticlesInEmCal();
  G4double totalEdepAllParticlesInHadCal = 
    theSteppingAction->getTotalEdepAllParticlesInHadCal();
  G4double totalEdepAllParticlesInMuonDetector = 
    theSteppingAction->getTotalEdepAllParticlesInMuonDetector();
  G4double exitingRadiusPrimaryMuon = 
    theSteppingAction->getExitingRadiusPrimaryMuon();

  int incidentParticleId = theSteppingAction->getPrimaryParticleId();
  G4double incidentParticleEnergy = theSteppingAction->getPrimaryParticleEnergy();

  if ( evt->GetEventID() == 0 ) {
    G4int isAlwaysKillLeadingHadronSet = 0;
    if ( getenv( "AlwaysKillLeadingHadron" ) ) {
      isAlwaysKillLeadingHadronSet = 1;
    }
    G4cout << "\t ----------------------------------" << G4endl
           << "\t Beam Particle PDG Id = " << incidentParticleId << G4endl
           << "\t Beam Particle Energy = " << incidentParticleEnergy / GeV
           << " GeV" << G4endl
           << "\t Is AlwaysKillLeadingHadron set? " 
           << isAlwaysKillLeadingHadronSet << G4endl
           << "\t ----------------------------------" << G4endl;
  }

  theSteppingAction->reset();

  G4cout << " ---  DetectorSliceEventAction::EndOfEventAction  ---  event= " 
	 << evt->GetEventID() << "   t=";
  eventTimer->Stop();
  G4double timeEventInSeconds = 0.0;
  if ( eventTimer->IsValid() ) {
    timeEventInSeconds = eventTimer->GetUserElapsed();
    G4cout << timeEventInSeconds << "s" << G4endl;
  } else {
    G4cout << "UNDEFINED" << G4endl;
  }

  // Fill the histograms/ntuple.
  DetectorSliceAnalysis* analysis = DetectorSliceAnalysis::getInstance();
  if ( analysis ) {
    analysis->fillNtuple( static_cast< float >( incidentParticleId ),
			  incidentParticleEnergy,
			  totalEdepAllParticlesInTracker, 
			  energyDepositedInActiveEmCalorimeterLayers,
			  totalEdepAllParticlesInEmCal, 
			  energyDepositedInActiveHadCalorimeterLayers,
			  totalEdepAllParticlesInHadCal, 
			  totalEdepAllParticlesInMuonDetector, 
                          exitingRadiusPrimaryMuon );
    analysis->endOfEvent( timeEventInSeconds );
  }

  // Now print out the information.
  //***LOOKHERE*** if ( evt->GetEventID() % 50 == 0 ) {  // Print info every 50 events.
  //if ( (evt->GetEventID() % 50) > -999 ) {
  //  G4cout << " DetectorSliceEventAction::EndOfEventAction : " << G4endl
  //	     << "\t Tracker : totalEdepAllParticlesInTracker              = " 
  //	     <<  totalEdepAllParticlesInTracker << G4endl
  //	     << "\t EM Cal :  totalEdepAllParticlesInEmCal                = "
  //         << totalEdepAllParticlesInEmCal << G4endl
  //         << "\t           energyDepositedInActiveEmCalorimeterLayers  = "
  //         << energyDepositedInActiveEmCalorimeterLayers << G4endl
  //	     << "\t HAD Cal :       totalEdepAllParticlesInHadCal         = "
  //         << totalEdepAllParticlesInHadCal << G4endl
  //         << "\t           energyDepositedInActiveHadCalorimeterLayers = "
  //         << energyDepositedInActiveHadCalorimeterLayers << G4endl
  //         << "\t Muon :    totalEdepAllParticlesInMuonDetector         = "
  //         << totalEdepAllParticlesInMuonDetector << G4endl 
  //         << "\t exitingRadiusPrimaryMuon                              = " 
  //         << exitingRadiusPrimaryMuon << G4endl;  //***DEBUG***
  //}

  // Extract the trajectories and draw them
  if ( G4VVisManager::GetConcreteInstance() ) {
    G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if ( trajectoryContainer ) n_trajectories = trajectoryContainer->entries();
    for ( G4int i=0; i<n_trajectories; i++ ) { 
      G4Trajectory* trj = ( G4Trajectory* )( (*( evt->GetTrajectoryContainer() ) )[i] );
      // trj->DrawTrajectory( 50 ); // Draw all tracks
      if ( trj->GetCharge() != 0. ) trj->DrawTrajectory( 50 ); // Draw only charged tracks
      // if ( trj->GetCharge() == 0. ) trj->DrawTrajectory( 50 ); // Draw only neutral tracks
      // if ( ( trj->GetCharge() == 0. ) && ( trj->GetParticleName() == "gamma" ) ) 
      //	trj->DrawTrajectory( 50 ); // Draw only gammas
      // if ( ( trj->GetCharge() == 0. ) && ( trj->GetParticleName() == "neutron" ) ) 
      //	trj->DrawTrajectory( 50 ); // Draw only neutrons
    }
  }

}


void DetectorSliceEventAction::instanciateSteppingAction() {      
  G4UserSteppingAction* theUserAction = const_cast< G4UserSteppingAction* >
    ( G4RunManager::GetRunManager()->GetUserSteppingAction() );
  if (theUserAction == 0) {
    theSteppingAction = new DetectorSliceSteppingAction;  
    G4RunManager::GetRunManager()->SetUserAction( theSteppingAction );
  }   
}

