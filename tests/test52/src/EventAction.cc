
#include "EventAction.hh"
#include "TargetLayerHitsCollection.hh"
#include "G4Event.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"


EventAction::EventAction() {

}


EventAction::~EventAction() {

}


void EventAction::BeginOfEventAction(const G4Event*) {

}


void EventAction::EndOfEventAction(const G4Event* event) {

  G4HCofThisEvent* hitsCollEvent = event -> GetHCofThisEvent();
  G4int nmbHitsCollections = 
           G4SDManager::GetSDMpointer() -> GetCollectionCapacity();

  for(G4int i=0; i < nmbHitsCollections; i++) {

     TargetLayerHitsCollection* hitsCollection =
  	        (TargetLayerHitsCollection*) hitsCollEvent -> GetHC(i);

     if(hitsCollection) {
    
        G4int nmbHits = hitsCollection -> entries();
        for(G4int j=0; j < nmbHits; j++) {
	
           Secondaries secParticles = 
                  ((*hitsCollection)[j]) -> GetSecParticles();

	   Secondaries::iterator iter = secParticles.begin();
	   Secondaries::iterator iter_end = secParticles.end();

           for(;iter != iter_end;iter++) {
              std::cout << iter -> first << " " << iter -> second << std::endl;
           }
        }
     }
  }

  if(G4VVisManager::GetConcreteInstance()) {

     G4HCofThisEvent* hitsCollEvent = event -> GetHCofThisEvent();
     G4int nmbHitsCollections = 
              G4SDManager::GetSDMpointer() -> GetCollectionCapacity();

     for(G4int i=0; i < nmbHitsCollections; i++) {

        TargetLayerHitsCollection* hitsCollection =
 	        (TargetLayerHitsCollection*) hitsCollEvent -> GetHC(i);

        if(hitsCollection) hitsCollection -> DrawAllHits();
     }

     G4TrajectoryContainer* trajContainer = event -> GetTrajectoryContainer();

     G4int nmbTrajectories = 0; 
     if (trajContainer) nmbTrajectories = trajContainer -> entries();
 
     for(G4int i=0; i < nmbTrajectories; i++) { 

        G4Trajectory* trajectory = 
           (G4Trajectory*) ((*(event->GetTrajectoryContainer()))[i]);
 
        if(trajectory) trajectory -> DrawTrajectory(50);
    }
  }

  G4int nmbEvents = event -> GetEventID();
  if(!(nmbEvents % 10000))
     std::cout << "INFORMATION: " << nmbEvents << " events processed." 
               << std::endl;

}
