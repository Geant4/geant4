#include "ML2EventAction.h"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"


CML2EventAction::CML2EventAction() :
  drawFlag("all" )
{
 }

 
CML2EventAction::~CML2EventAction()
{
 }
 
void CML2EventAction::BeginOfEventAction(const G4Event*)
{
}

 
void CML2EventAction::EndOfEventAction(const G4Event* evt)
{  
 // extract the trajectories and draw them ...

  if (G4VVisManager::GetConcreteInstance())
    {
      G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
      G4int n_trajectories = 0;
      if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

      for (G4int i=0; i<n_trajectories; i++) 
        {
			G4Trajectory* trj = (G4Trajectory*)
			((*(evt->GetTrajectoryContainer()))[i]);
			if(drawFlag == "all") trj->DrawTrajectory(50);
			else if((drawFlag == "charged")&&(trj->GetCharge() != 0.))
			trj->DrawTrajectory(50);
			else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.))
			trj->DrawTrajectory(50);	


		}
    }
 }
