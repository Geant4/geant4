//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
#include "LXeEventAction.hh"
#include "LXeScintHit.hh"
#include "LXePMTHit.hh"
#include "LXeUserEventInformation.hh"
#include "LXeTrajectory.hh"
#include "RecorderBase.hh"

#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UImanager.hh"

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXeEventAction::LXeEventAction(RecorderBase* r)
  :recorder(r),saveThreshold(0),scintCollID(-1),pmtCollID(-1),verbose(0),
   pmtThreshold(1),forcedrawphotons(false),forcenophotons(false)
{
  eventMessenger=new LXeEventMessenger(this);
}
 
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
LXeEventAction::~LXeEventAction(){}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeEventAction::BeginOfEventAction(const G4Event* anEvent){
  
  //New event, add the user information object
  G4EventManager::
    GetEventManager()->SetUserInformation(new LXeUserEventInformation);
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if(scintCollID<0)
    scintCollID=SDman->GetCollectionID("scintCollection");
  if(pmtCollID<0)
    pmtCollID=SDman->GetCollectionID("pmtHitCollection");

  if(recorder)recorder->RecordBeginOfEvent(anEvent);
}
 
//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeEventAction::EndOfEventAction(const G4Event* anEvent){
  
  LXeUserEventInformation* eventInformation
    =(LXeUserEventInformation*)anEvent->GetUserInformation();
  
  G4TrajectoryContainer* trajectoryContainer=anEvent->GetTrajectoryContainer();
  
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  // extract the trajectories and draw them
  if (G4VVisManager::GetConcreteInstance()){
    for (G4int i=0; i<n_trajectories; i++){ 
      LXeTrajectory* trj = (LXeTrajectory*)
	((*(anEvent->GetTrajectoryContainer()))[i]);
      if(trj->GetParticleName()=="opticalphoton"){
	trj->SetForceDrawTrajectory(forcedrawphotons);
	trj->SetForceNoDrawTrajectory(forcenophotons);
      }
      trj->DrawTrajectory(50);
    }
  }
  
  LXeScintHitsCollection* SHC = 0;
  LXePMTHitsCollection* PHC = 0;
  G4HCofThisEvent* HCE = anEvent->GetHCofThisEvent();
  
  //Get the hit collections
  if(HCE){
    if(scintCollID>=0)SHC = (LXeScintHitsCollection*)(HCE->GetHC(scintCollID));
    if(pmtCollID>=0)PHC = (LXePMTHitsCollection*)(HCE->GetHC(pmtCollID));
  }

  //Hits in scintillator
  if(SHC){
    int n_hit = SHC->entries();
    G4ThreeVector  eWeightPos(0.);     
    G4double edep;
    G4double edepMax=0;

    for(int i=0;i<n_hit;i++){ //gather info on hits in scintillator
      edep=(*SHC)[i]->GetEdep();
      eventInformation->IncEDep(edep); //sum up the edep
      eWeightPos += (*SHC)[i]->GetPos()*edep;//calculate energy weighted pos
      if(edep>edepMax){
	edepMax=edep;//store max energy deposit
	G4ThreeVector posMax=(*SHC)[i]->GetPos();
	eventInformation->SetPosMax(posMax,edep);
      }
    }
    if(eventInformation->GetEDep()==0.){
      if(verbose>0)G4cout<<"No hits in the scintillator this event."<<G4endl;
    }    
    else{
      //Finish calculation of energy weighted position
      eWeightPos/=eventInformation->GetEDep();
      eventInformation->SetEWeightPos(eWeightPos);
      if(verbose>0){
	G4cout << "\tEnergy weighted position of hits in LXe : "
	       << eWeightPos/mm << G4endl;
      }
    }
    if(verbose>0){
    G4cout << "\tTotal energy deposition in scintillator : "
	   << eventInformation->GetEDep() / keV << " (keV)" << G4endl;
    }
  }
  
  if(PHC){
    G4ThreeVector reconPos(0.,0.,0.);
    G4int pmts=PHC->entries();
    //Gather info from all PMTs
    for(G4int i=0;i<pmts;i++){
      eventInformation->IncHitCount((*PHC)[i]->GetPhotonCount());
      reconPos+=(*PHC)[i]->GetPMTPos()*(*PHC)[i]->GetPhotonCount();
      if((*PHC)[i]->GetPhotonCount()>=pmtThreshold){
	eventInformation->IncPMTSAboveThreshold();
      }
      else{//wasnt above the threshold, turn it back off
	(*PHC)[i]->SetDrawit(false);
      }
    }
    
    if(eventInformation->GetHitCount()>0){//dont bother unless there were hits
      reconPos/=eventInformation->GetHitCount();
      if(verbose>0){
	G4cout << "\tReconstructed position of hits in LXe : "
	       << reconPos/mm << G4endl;
      }
      eventInformation->SetReconPos(reconPos);
    }
    PHC->DrawAllHits();
  }

  if(verbose>0){
    //End of event output. later to be controlled by a verbose level
    G4cout << "\tNumber of photons that hit PMTs in this event : "
	   << eventInformation->GetHitCount() << G4endl;
    G4cout << "\tNumber of PMTs above threshold("<<pmtThreshold<<") : "
	   << eventInformation->GetPMTSAboveThreshold() << G4endl;
    G4cout << "\tNumber of photons produced by scintillation in this event : "
	   << eventInformation->GetPhotonCount_Scint() << G4endl;
    G4cout << "\tNumber of photons produced by cerenkov in this event : "
	   << eventInformation->GetPhotonCount_Ceren() << G4endl;
    G4cout << "\tNumber of photons absorbed (OpAbsorption) in this event : "
	   << eventInformation->GetAbsorptionCount() << G4endl;
    G4cout << "\tNumber of photons absorbed at boundaries (OpBoundary) in "
	   << "this event : " << eventInformation->GetBoundaryAbsorptionCount() 
	   << G4endl;
    G4cout << "Unacounted for photons in this event : " 
	   << (eventInformation->GetPhotonCount_Scint() + 
	       eventInformation->GetPhotonCount_Ceren() - 
	       eventInformation->GetAbsorptionCount() -
	       eventInformation->GetHitCount() - 
	       eventInformation->GetBoundaryAbsorptionCount())
	   << G4endl;
  }
  //If we have set the flag to save 'special' events, save here
  if(saveThreshold&&eventInformation->GetPhotonCount() <= saveThreshold)
    G4RunManager::GetRunManager()->rndmSaveThisEvent();

  if(recorder)recorder->RecordEndOfEvent(anEvent);
}

//_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
void LXeEventAction::SetSaveThreshold(G4int save){
/*Sets the save threshold for the random number seed. If the number of photons
generated in an event is lower than this, then save the seed for this event
in a file called run###evt###.rndm
*/
  saveThreshold=save;
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4RunManager::GetRunManager()->SetRandomNumberStoreDir("random/");
  //  G4UImanager::GetUIpointer()->ApplyCommand("/random/setSavingFlag 1");
}




