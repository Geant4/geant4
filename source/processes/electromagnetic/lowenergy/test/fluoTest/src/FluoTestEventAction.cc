//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FluoTestEventAction.hh"
#include "FluoTestSensorHit.hh"
#include "FluoTestEventActionMessenger.hh"

#ifdef G4ANALYSIS_USE
#include "FluoTestAnalysisManager.hh"
#endif

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
#include "G4UnitsTable.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4ANALYSIS_USE
FluoTestEventAction::FluoTestEventAction(FluoTestAnalysisManager* aMgr):
  drawFlag("all"),
  HPGeCollID(-1),
  eventMessenger(0), 
  printModulo(1),fAnalysisManager(aMgr)
 {
   eventMessenger = new FluoTestEventActionMessenger(this);

}

#else

FluoTestEventAction::FluoTestEventAction()
  :drawFlag("all"),
   HPGeCollID(-1),
   eventMessenger(0),
   printModulo(1)
 
{
  eventMessenger = new FluoTestEventActionMessenger(this);
}

#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


FluoTestEventAction::~FluoTestEventAction()
{
   delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestEventAction::BeginOfEventAction(const G4Event* evt)
{
  
  /*G4int evtNb = evt->GetEventID(); // Returns the event ID
     if (evtNb%printModulo == 0)
   { 
  G4cout << "\n---> Begin of event: " << evtNb << G4endl;
  HepRandom::showEngineStatus();
      }
  */
  if (HPGeCollID==-1)
    {
      G4SDManager * SDman = G4SDManager::GetSDMpointer();
      HPGeCollID = SDman->GetCollectionID("HPGeCollection");
      //the pointer points to the ID number of the sensitive detector
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestEventAction::EndOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID(); 
  
 // extracted from hits, compute the total energy deposit (and total charged
  // track length) 
  
  
    G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
    
    FluoTestSensorHitsCollection* HPGeHC = 0;
    G4int n_hit = 0;
    G4double totEnergyDetect=0., totEnergy=0, energyD=0.;
    
    if (HCE) HPGeHC = (FluoTestSensorHitsCollection*)(HCE->GetHC(HPGeCollID));
    if(HPGeHC)
      {
	n_hit = HPGeHC->entries();
	for (G4int i=0;i<n_hit;i++)
	  {
	    totEnergy += (*HPGeHC)[i]->GetEdepTot(); 
	    /*
#ifdef G4ANALYSIS_USE
	    fAnalysisManager->InsGamDet((*HPGeHC)[i]->GetEdepTot()/keV);  
#endif;
	    */ 
	    energyD    = (*HPGeHC)[i]->RandomCut();
	    totEnergyDetect += energyD;
	 
#ifdef G4ANALYSIS_USE
	    fAnalysisManager->InsGamDet(energyD/keV);  
#endif
	  }
      }
    // print this information event by event (modulo n)  	
    
    if (evtNb%printModulo == 0) 
      { 
    if (totEnergy!=0)
      {
	// G4cout << "---> End of event: " << evtNb << G4endl;
	//G4cout << " total energy detected: " << G4std::setw(14) 
	// << G4BestUnit(totEnergyDetect,"Energy");
    //G4cout << "\n     " << n_hit
	//      << " hits are stored in HPGeCollection." << G4endl;
	/*
#ifdef G4ANALYSIS_USE
	fAnalysisManager->InsDetETot(totEnergyDetect/keV);  
#endif; 
	*/
      }
    
      }
    /*   
#ifdef G4ANALYSIS_USE
    fAnalysisManager->EndOfEvent(evtNb);
#endif 
    */ 
  // extract the trajectories and draw them
    
  if (G4VVisManager::GetConcreteInstance())
    {
      
      G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
      G4int n_trajectories = 0;
      if (trajectoryContainer) n_trajectories = trajectoryContainer->size();
      
      for (G4int i=0; i<n_trajectories; i++) 
	{ G4Trajectory* trj = (G4Trajectory*)((*(evt->GetTrajectoryContainer()))[i]);
	if (drawFlag == "all") trj->DrawTrajectory(50);
	else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
	  trj->DrawTrajectory(50);
	else if ((drawFlag == "neutral")&&(trj->GetCharge() == 0.))
	  trj->DrawTrajectory(50);				   
	}
    }             
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




