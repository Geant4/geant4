//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// Rich advanced example for Geant4
// RichTbEventAction.cc for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#include "RichTbEventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4VVisManager.hh"
#include "RichTbVisManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "RichTbSD.hh"
#include "RichTbHit.hh"
#include "RichTbRunConfig.hh"
#include "RichTbAnalysisManager.hh"
#include "RichTbIOData.hh"

RichTbEventAction::RichTbEventAction(){
  RichTbCollID = -1;

}

RichTbEventAction::RichTbEventAction(RichTbRunConfig* RConfig, 
    RichTbVisManager* RVisManager, RichTbIOData* RIOData){
  RichTbCollID = -1;
  runConfiguration=RConfig;

  pVisManager=RVisManager;
  rTbIOData = RIOData;
  
  
}
RichTbEventAction::~RichTbEventAction(){


}
void RichTbEventAction::BeginOfEventAction(const G4Event* evt){




  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  if(RichTbCollID<0){



    RichTbCollID = SDman->GetCollectionID("RichTbHitsCollection");

  }

#ifdef G4ANALYSIS_USE
  RichTbAnalysisManager * analysis = RichTbAnalysisManager::getInstance();
  analysis->BeginOfEventAnalysis(evt);

#endif


}
void RichTbEventAction::EndOfEventAction(const G4Event* evt){

  RichTbRunConfig* RConfig =  runConfiguration;

  // first get the trajectories
  G4TrajectoryContainer* trajectoryContainer=evt->GetTrajectoryContainer();
  G4int n_trajectories=0;

  G4int totalhits = 0;

  if(trajectoryContainer){n_trajectories=trajectoryContainer->entries();
    G4cout << "     " << n_trajectories
         << " Tracks are stored in Trajectorycontainer." << G4endl;
  }
    // Now get the hits

  if(RichTbCollID<0) return;


  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  RichTbHitsCollection* RHC = 0;
  G4int n_hit = 0; 

  if(HCE)
  {

    RHC = (RichTbHitsCollection*)(HCE->GetHC(RichTbCollID));

  }

  

  if(RHC)
  {



    n_hit = RHC->entries();


    for(G4int i=0;i<n_hit;i++)

      {

	totalhits += i;

      }




    //    G4cout << "     " << n_hit
    //   << " hits are stored in RichTbHitsCollection." << G4endl;
  }

    if(RConfig->getRichTbDrawTrajectory() > 0){    
     for (G4int i=0; i<n_trajectories; i++){
     (*(evt->GetTrajectoryContainer()))[i]->DrawTrajectory();
     }
   }
  

    //now for the hits.

    if(RHC) {
     G4int n_hith = RHC->entries();
    



     for (G4int ih=0; ih<n_hith; ih++ ) {
       RichTbHit* aHit = (*RHC)[ih];
       aHit->DrawWithVisM(pVisManager);
     }

    }
  //Now to Fill the Histograms

#ifdef G4ANALYSIS_USE
    RichTbAnalysisManager * analysis = RichTbAnalysisManager::getInstance();
    analysis->EndOfEventAnalysis(evt);
#endif

   // Now to write out the data
   if(RConfig -> DoWriteOutputFile() ) {

   rTbIOData -> WriteOutEventHeaderData(evt);
   rTbIOData -> WriteOutHitData(evt);
   G4cout << ">>> End of Event  Number  " << evt->GetEventID() << G4endl;

   }

}









