#include "MyEventAction.hh"

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
#include "MySensitiveCalorimeter.hh"
#include "MyCalorimeterHit.hh"

#include "MySteppingAction.hh"
#include "G4RunManager.hh"

#include "MyAnalysis.hh"


MyEventAction::MyEventAction() : theSteppingAction(0) {
  instanciateSteppingAction();
}


MyEventAction::~MyEventAction() {}


void MyEventAction::StartOfEventAction(const G4Event* evt) { 
  G4cout << "\n---> Begin of event: " << evt->GetEventID() << G4endl;
}


void MyEventAction::EndOfEventAction(const G4Event* evt){

  // Get the IDs of the collections of calorimeter hits and tracker hits. 
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  G4int calorimeterCollID = SDman->GetCollectionID(colNam="calCollection");

  // Extract the hits: first for the calorimeter.
  G4HCofThisEvent* collectionEventHits = evt->GetHCofThisEvent();
  MyCalorimeterHitsCollection* collectionCaloHits = 0;
  G4int numCalorimeterHits = 0;
  G4double energyDepositedInActiveCalorimeterLayer = 0.0;
  if ( collectionEventHits ) {
    collectionCaloHits = dynamic_cast< MyCalorimeterHitsCollection* >
      ( collectionEventHits->GetHC(calorimeterCollID) );
  }
  if ( collectionCaloHits ) {
    numCalorimeterHits = collectionCaloHits->entries();
    for ( G4int i=0; i < numCalorimeterHits; i++ ) {
      energyDepositedInActiveCalorimeterLayer += (*collectionCaloHits)[i]->GetEdep(); 
    }
  }

  // Get the information about the total release of energy in all the calorimeter
  // (not only in the sensitive parts) and the Id and kinetic energy of the 
  // incident (primary) particle from the Stepping Action.
  G4double totalEdepAllParticles = theSteppingAction->getTotalEdepAllParticles();
  int incidentParticleId = theSteppingAction->getPrimaryParticleId();
  G4double incidentParticleEnergy = theSteppingAction->getPrimaryParticleEnergy();

  theSteppingAction->reset();

  // Fill the histograms/ntuple.
  MyAnalysis* analysis = MyAnalysis::getInstance();
  analysis->fillNtuple( static_cast< float >( incidentParticleId ),
			incidentParticleEnergy / MeV,
			energyDepositedInActiveCalorimeterLayer / MeV,
			totalEdepAllParticles / MeV );

  // Now print out the information.
  G4cout << " ---  MyEventAction::EndOfEventAction  ---    event = " 
	 << evt->GetEventID() << G4endl;
  if ( evt->GetEventID() % 50 == 0 ) {  // Print info every 50 events.
    if ( numCalorimeterHits ) {
      G4cout << "\t Energy deposited in the Active calorimeter layer = " 
             << energyDepositedInActiveCalorimeterLayer / MeV << " (MeV)" 
             << G4endl;  //***DEBUG***
    } else {
      G4cout << "\t NO HIT in the Active calorimeter layer." << G4endl; //***DEBUG***
    }
    G4cout << "\t totalEdepAllParticles     = " << totalEdepAllParticles / GeV 
	   << " (GeV)" << G4endl;  //***DEBUG***
  }  

  // Extract the trajectories and draw them
  if ( G4VVisManager::GetConcreteInstance() ) {
    G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if ( trajectoryContainer ) n_trajectories = trajectoryContainer->entries();
    for ( G4int i=0; i<n_trajectories; i++ ) { 
      G4Trajectory* trj = (G4Trajectory*) ((*(evt->GetTrajectoryContainer()))[i]);
      // trj->DrawTrajectory(50); // Draw all tracks
      if ( trj->GetCharge() != 0. ) trj->DrawTrajectory(50); // Draw only charged tracks
      // if ( trj->GetCharge() == 0. ) trj->DrawTrajectory(50); // Draw only neutral tracks
      // if ( ( trj->GetCharge() == 0. ) && ( trj->GetParticleName() == "gamma" ) ) 
      //	trj->DrawTrajectory(50); // Draw only gammas
      // if ( ( trj->GetCharge() == 0. ) && ( trj->GetParticleName() == "neutron" ) ) 
      //	trj->DrawTrajectory(50); // Draw only neutrons
    }
  }

}


void MyEventAction::instanciateSteppingAction() {      
  G4UserSteppingAction* theUserAction = const_cast<G4UserSteppingAction*>
    ( G4RunManager::GetRunManager()->GetUserSteppingAction() );
  if (theUserAction == 0) {
    theSteppingAction = new MySteppingAction;  
    G4RunManager::GetRunManager()->SetUserAction(theSteppingAction);
  }   
}

