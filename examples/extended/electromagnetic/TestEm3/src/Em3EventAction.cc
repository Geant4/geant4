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
// $Id: Em3EventAction.cc,v 1.12 2001-11-28 17:54:46 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em3EventAction.hh"

#include "Em3RunAction.hh"
#include "Em3PrimaryGeneratorAction.hh"
#include "Em3DetectorConstruction.hh"
#include "Em3CalorHit.hh"
#include "Em3EventActionMessenger.hh"
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

#ifndef G4NOHIST
 #include "CLHEP/Hist/HBookFile.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em3EventAction::Em3EventAction(Em3RunAction* run,Em3PrimaryGeneratorAction* kin,
                               Em3DetectorConstruction* det)
:Em3Run(run),Em3Kin(kin),Detector(det),calorimeterCollID(-1),drawFlag("none"),
 printModulo(10000),eventMessenger(0)
{
  eventMessenger = new Em3EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em3EventAction::~Em3EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em3EventAction::BeginOfEventAction(const G4Event* evt)
{   
 G4int evtNb = evt->GetEventID();
 
 //survey printing
 if (evtNb%printModulo == 0) 
    G4cout << "\n---> Begin Of Event: " << evtNb << G4endl;

 // initialize Hits collection    
 if (calorimeterCollID==-1)
  {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    calorimeterCollID = SDman->GetCollectionID("CalCollection");
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em3EventAction::EndOfEventAction(const G4Event* evt)
{ 
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  Em3CalorHitsCollection* CHC = NULL;
  G4int NbHits=0;
  G4int NbOfAbsor=Detector->GetNbOfAbsor();
  G4double Ebeam = Em3Kin->GetParticleGun()->GetParticleEnergy();
  G4double totEAbs, totLAbs;
       
  if (HCE) CHC = (Em3CalorHitsCollection*)(HCE->GetHC(calorimeterCollID));

  if (CHC)
    {
     NbHits = CHC->entries();
     for (G4int k=0; k<NbOfAbsor; k++)
        {totEAbs=totLAbs=0.;
         for (G4int j=0;j<NbHits;j++)
            {
	     totEAbs += (*CHC)[j]->GetEdepAbs(k); 
             totLAbs += (*CHC)[j]->GetTrakAbs(k);     
            }
         Em3Run->fillPerEvent(k,totEAbs,totLAbs);
	 
#ifndef G4NOHIST        
         //fill histo
         //	 
	 Em3Run->GetHisto(k)->accumulate(totEAbs/Ebeam);
#endif  	 	 
       }
    }
    
  if (G4VVisManager::GetConcreteInstance())
    {
     G4TrajectoryContainer* trajectoryContainer = evt->GetTrajectoryContainer();
     G4int n_trajectories = 0;
     if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();

     for (G4int i=0; i<n_trajectories; i++) 
        { G4Trajectory* trj = (G4Trajectory*)
	                             ((*(evt->GetTrajectoryContainer()))[i]);
          if (drawFlag == "all") trj->DrawTrajectory(50);
          else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
                                  trj->DrawTrajectory(50); 
        }
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


