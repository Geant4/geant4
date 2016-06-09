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
//
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//    ****************************************************
//    *      UltraEventAction.cc
//    ****************************************************
//
//    Ultra EventAction class. The UltraAnalysisManager class is used for histogram
//    filling if the G4ANALYSIS_USE environment variable is set.
//
#include "UltraEventAction.hh"
#include "UltraRunAction.hh"
#include "UltraPrimaryGeneratorAction.hh"
#include "UltraOpticalHit.hh"

#include "G4RunManager.hh" 
#include "G4Event.hh"
#include "G4EventManager.hh" 
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4GeneralParticleSource.hh" 

#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"

#ifdef G4ANALYSIS_USE
#include "UltraAnalysisManager.hh"
#endif
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraEventAction::UltraEventAction(UltraRunAction* run)
  :UltraRun(run),OpticalHitsCollID(-1)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

UltraEventAction::~UltraEventAction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraEventAction::BeginOfEventAction(const G4Event* evt)
{
  G4int printModulo = 100;

  evtNb = evt->GetEventID();

  G4SDManager * SDman = G4SDManager::GetSDMpointer(); 


  if(OpticalHitsCollID==-1) {
    OpticalHitsCollID = SDman->GetCollectionID("OpticalHitsCollection");
  }


  if (evtNb%printModulo == 0)
    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void UltraEventAction::EndOfEventAction(const G4Event* evt)
{

G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
UltraOpticalHitsCollection* OpticalHitsColl = 0;
  
if(HCE){

  if(OpticalHitsCollID != -1) OpticalHitsColl = 
  (UltraOpticalHitsCollection*)(HCE->GetHC(OpticalHitsCollID));

}
  G4int nOptHits = 0 ; 

  if(OpticalHitsColl){

    nOptHits = OpticalHitsColl->entries();

#ifdef ULTRA_VERBOSE
    if (nOptHits > 0){
     G4cout << " Optical Hit # " << " " << "Energy (eV)" <<  " " << "x,y,z (cm)" << G4endl ;
    }
#endif

    for(G4int iHit=0; iHit<nOptHits; iHit++){

      G4double HitEnergy ;
      HitEnergy = (*OpticalHitsColl)[iHit]->GetEnergy() ;
      G4ThreeVector HitPos = (*OpticalHitsColl)[iHit]->GetPosition()/cm;

#ifdef G4ANALYSIS_USE
      UltraAnalysisManager* analysis = UltraAnalysisManager::getInstance();
      analysis->FillHistogram(1,HitEnergy/eV);
#endif

    }

  }

#ifdef G4ANALYSIS_USE
  UltraAnalysisManager* analysis = UltraAnalysisManager::getInstance();
  analysis->FillHistogram(2,nOptHits);
#endif

  // Define visualization atributes:

  G4String drawFlag = "all";
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();

  if(pVVisManager)
  {
   G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
   G4int n_trajectories = 0;
   if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();  
   for(G4int i=0; i<n_trajectories; i++) 
      { G4Trajectory* trj = (G4Trajectory *)((*(evt->GetTrajectoryContainer()))[i]);
        if (drawFlag == "all") trj->DrawTrajectory();
        else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                               trj->DrawTrajectory(); 
      }
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
