// -------------------------------------------------------------------
// $Id: MicrobeamEventAction.cc,v 1.2 2006-04-10 14:47:32 sincerti Exp $
// -------------------------------------------------------------------

#include "G4Event.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "Randomize.hh"

#include "MicrobeamEventAction.hh"
#include "MicrobeamRunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamEventAction::MicrobeamEventAction(MicrobeamRunAction* run)
  :Run(run),drawFlag("all"),printModulo(10000)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicrobeamEventAction::~MicrobeamEventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrobeamEventAction::BeginOfEventAction(const G4Event* evt)
{  
  G4int evtNb = evt->GetEventID();
  Run->SetNumEvent(evtNb);
  Run->SetDoseN(0);
  Run->SetDoseC(0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicrobeamEventAction::EndOfEventAction(const G4Event* evt)
{  
 
// SAVE TOTAL ABSORBED DOSE IN PHANTOM

if (Run->GetDoseN()>0 || Run->GetDoseC()>0) 
{
	FILE* myFile;
	myFile=fopen("dose.txt","a");
	fprintf(myFile,"%e %e\n",Run->GetDoseN(),Run->GetDoseC());
	fclose (myFile);
        G4cout << "   ===> The incident alpha particle has reached the targeted cell :" << G4endl;
	G4cout << "   -----> total absorbed dose within Nucleus   is (Gy) = " << Run->GetDoseN() << G4endl;
	G4cout << "   -----> total absorbed dose within Cytoplasm is (Gy) = " << Run->GetDoseC() << G4endl;
	G4cout << G4endl;
}
else
{
	G4cout << "   ===> Sorry, the incident alpha particle has missed the targeted cell !" << G4endl;
	G4cout << G4endl;
}


if (G4VVisManager::GetConcreteInstance())
    {
      G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
      G4int n_trajectories = 0;
      if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
      
      for (G4int i=0; i<n_trajectories; i++) 
        { 
	  //G4cout<< "Iteration" << i <<G4endl;
	  G4Trajectory* trj = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
	  if (drawFlag == "all")
	    { 
	      trj->DrawTrajectory(50);
	    }
	  else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))	    
	    {
	      trj->DrawTrajectory(50);
	    }
	} 

            //save rndm status
      if (Run->GetRndmFreq() == 2)
	{
	  HepRandom::saveEngineStatus("endOfEvent.rndm");   
	  G4int evtNb = evt->GetEventID();
	  if (evtNb%printModulo == 0)
	    { 
	      //G4cout << "\n---> End of Event: " << evtNb << G4endl;
	      //HepRandom::showEngineStatus();
	    }
	}    
    }

}
