#include "RandomCaloEventAction.hh"

#include "G4Event.hh"
#include "G4ios.hh"

#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"

#include "RandomCaloSteppingAction.hh"
#include "G4RunManager.hh"


RandomCaloEventAction::RandomCaloEventAction() : 
  theSteppingAction( 0 ),
  numberOfEvents( 0 ),
  sumTotalDepositedEnergy( 0.0 ) {
  instanciateSteppingAction();
}


RandomCaloEventAction::~RandomCaloEventAction() {
  if ( numberOfEvents > 0 ) {
    G4cout << G4endl
           << " ----------------------------------- " << G4endl
           << " Average deposited energy = " 
	   << sumTotalDepositedEnergy / static_cast< G4double >( numberOfEvents )
           << "  MeV " << G4endl
           << " ----------------------------------- " << G4endl
           << G4endl;
  }
}


void RandomCaloEventAction::StartOfEventAction( const G4Event* ) { 
  //G4cout << "\n---> Begin of event: " << evt->GetEventID() << G4endl;
}


void RandomCaloEventAction::EndOfEventAction( const G4Event* evt ) {

  G4cout << " ---  RandomCaloEventAction::EndOfEventAction  ---    event = " 
	 << evt->GetEventID() << G4endl;

  G4double totalEdepAllParticles = theSteppingAction->getTotalEdepAllParticles();

  sumTotalDepositedEnergy += totalEdepAllParticles;

  numberOfEvents++;

  theSteppingAction->reset();

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


void RandomCaloEventAction::instanciateSteppingAction() {      
  G4UserSteppingAction* theUserAction = const_cast< G4UserSteppingAction* >
    ( G4RunManager::GetRunManager()->GetUserSteppingAction() );
  if (theUserAction == 0) {
    theSteppingAction = new RandomCaloSteppingAction;  
    G4RunManager::GetRunManager()->SetUserAction( theSteppingAction );
  }   
}

