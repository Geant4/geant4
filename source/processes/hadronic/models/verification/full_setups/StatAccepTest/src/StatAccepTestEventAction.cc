#include "StatAccepTestEventAction.hh"

#include "G4Event.hh"
#include "G4ios.hh"

#ifdef G4VIS_USE
 #include "G4TrajectoryContainer.hh"
 #include "G4Trajectory.hh"
 #include "G4VVisManager.hh"
#endif

#include "G4SDManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "StatAccepTestSensitiveCalorimeter.hh"
#include "StatAccepTestCalorimeterHit.hh"

#include "StatAccepTestSteppingAction.hh"
#include "G4RunManager.hh"
#include "G4Timer.hh"

#include "StatAccepTestAnalysis.hh"


StatAccepTestEventAction::StatAccepTestEventAction() : 
  theSteppingAction( 0 ) {
  instanciateSteppingAction();
  eventTimer = new G4Timer;
}


StatAccepTestEventAction::~StatAccepTestEventAction() {
  delete eventTimer;
}


void StatAccepTestEventAction::BeginOfEventAction( const G4Event* ) {
  //G4cout << "\n---> Begin of event: " << evt->GetEventID() << G4endl;
  eventTimer->Start();
}


void StatAccepTestEventAction::EndOfEventAction( const G4Event* evt ) {

  // Get the IDs of the collections of calorimeter hits and tracker hits. 
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  G4int calorimeterCollID = SDman->GetCollectionID( colNam="calCollection" );

  // Extract the hits: first for the calorimeter.
  G4HCofThisEvent* collectionEventHits = evt->GetHCofThisEvent();
  StatAccepTestCalorimeterHitsCollection* collectionCaloHits = 0;
  G4int numCalorimeterHits = 0;
  G4double energyDepositedInActiveCalorimeterLayer = 0.0;
  if ( collectionEventHits ) {
    collectionCaloHits = dynamic_cast< StatAccepTestCalorimeterHitsCollection* >
      ( collectionEventHits->GetHC( calorimeterCollID ) );
  }
  if ( collectionCaloHits ) {
    numCalorimeterHits = collectionCaloHits->entries();
    for ( G4int i=0; i < numCalorimeterHits; i++ ) {
      energyDepositedInActiveCalorimeterLayer += (*collectionCaloHits)[i]->GetEdep(); 
      //G4cout << " StatAccepTestEventAction::EndOfEventAction : " << G4endl
      //       << " \t iHit=" << i 
      //       << "  layer="  << (*collectionCaloHits)[i]->GetLayer()
      //       << "  energy=" << (*collectionCaloHits)[i]->GetEdep()
      //       << G4endl;           //***DEBUG***
    }
  }

  //G4cout << " StatAccepTestEventAction::EndOfEventAction : numCalorimeterHits = " 
  //       << numCalorimeterHits << std::endl; //***DEBUG***

  // Get the information about the total release of energy in all the calorimeter
  // (not only in the sensitive parts) and the Id and kinetic energy of the 
  // incident (primary) particle from the Stepping Action.
  G4double totalEdepAllParticles = theSteppingAction->getTotalEdepAllParticles();
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

  //G4cout << " ---  StatAccepTestEventAction::EndOfEventAction  ---  event= " 
  G4cout << " ---  EndOfEventAction  ---  event= " 
	 << evt->GetEventID() << "   t=";
  eventTimer->Stop();
  G4double timeEventInSeconds = 0.0;
  if ( eventTimer->IsValid() ) {
    timeEventInSeconds = eventTimer->GetUserElapsed();
    G4cout << timeEventInSeconds << "s";
  } else {
    G4cout << "UNDEFINED";
  }

  //G4cout << "   evis=" << energyDepositedInActiveCalorimeterLayer;

  G4cout << "   random=" << CLHEP::HepRandom::getTheEngine()->flat() << G4endl;

  // Fill the histograms/ntuple.
  StatAccepTestAnalysis* analysis = StatAccepTestAnalysis::getInstance();
  if ( analysis ) {
    analysis->fillNtuple( static_cast< float >( incidentParticleId ),
			  incidentParticleEnergy / MeV,
			  energyDepositedInActiveCalorimeterLayer / MeV,
			  totalEdepAllParticles / MeV );
    analysis->endOfEvent( timeEventInSeconds );
  }

  // Now print out the information.
  //***LOOKHERE*** if ( evt->GetEventID() % 50 == 0 ) {  // Print info every 50 events.
  if ( (evt->GetEventID() % 50) > -999 ) {
    if ( numCalorimeterHits ) {
      //G4cout << " StatAccepTestEventAction::EndOfEventAction : " << G4endl
      //       << "\t Energy deposited in the Active calorimeter layer = " 
      //       << energyDepositedInActiveCalorimeterLayer / MeV << " (MeV)" 
      //       << G4endl;  //***DEBUG***
    } else {
      //G4cout << " StatAccepTestEventAction::EndOfEventAction : " << G4endl
      //       << "\t NO HIT in the Active calorimeter layer." << G4endl; //***DEBUG***
    }
    //G4cout << "\t totalEdepAllParticles     = " << totalEdepAllParticles / MeV 
    //       << " MeV" << G4endl;  //***DEBUG***
  }  

#ifdef G4VIS_USE
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
#endif

}


void StatAccepTestEventAction::instanciateSteppingAction() {      
  G4UserSteppingAction* theUserAction = const_cast< G4UserSteppingAction* >
    ( G4RunManager::GetRunManager()->GetUserSteppingAction() );
  if (theUserAction == 0) {
    theSteppingAction = new StatAccepTestSteppingAction;  
    G4RunManager::GetRunManager()->SetUserAction( theSteppingAction );
  }   
}

