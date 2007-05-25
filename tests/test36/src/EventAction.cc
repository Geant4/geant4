#include "EventAction.hh"

#include "DetectorConstruction.hh"
#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(DetectorConstruction* det)
:detector(det)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event* )
{   

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{


// print modulo
  G4int event_id = evt->GetEventID();
  if ( (event_id + 1 )%1000000 == 0) G4cout << (event_id + 1 )/1000000.<<"  million of events" << G4endl;

// visualization
  if (G4VVisManager::GetConcreteInstance())
    {
     G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
     G4int n_trajectories = 0;
     if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
     for (G4int i=0; i<n_trajectories; i++) 
        { G4Trajectory* trj = (G4Trajectory*)
	                             ((*(evt->GetTrajectoryContainer()))[i]);
          if (drawFlag == "all") trj->DrawTrajectory(1000);
          else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                                  trj->DrawTrajectory(1000); 
          else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.))
                                 trj->DrawTrajectory(1000); 
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


