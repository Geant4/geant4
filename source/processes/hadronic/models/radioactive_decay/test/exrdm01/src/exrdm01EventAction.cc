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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		SSAEventAction.cc
//
// Version:		0.b.4
// Date:		16/08/99
// Author:		F Lei
// Organisation:	DERA UK
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		12115/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// CHANGE HISTORY
// --------------
//
// 16 August 1999, F Lei, DERA UK
// Adapted from a verson by Bill Lockman, SLAC, to whom all credits go:
//
// $Id: exrdm01EventAction.cc,v 1.1 2003-10-08 16:31:47 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//#include <rw/tvordvec.h>
#include "G4ios.hh"
#include "exrdm01EventActionMessenger.hh"
#include "exrdm01EventAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"

extern G4bool drawEvent;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdm01EventAction::exrdm01EventAction()
  : drawFlag("all"),eventMessenger(NULL)
{
  eventMessenger = new exrdm01EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

exrdm01EventAction::~exrdm01EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdm01EventAction::BeginOfEventAction(const G4Event* Ev)
{drawEvent=false;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void exrdm01EventAction::EndOfEventAction(const G4Event* Ev)
{
  if (drawEvent){
    const G4Event* evt = fpEventManager->GetConstCurrentEvent();

    //  G4cout << ">>> Event " << evt->GetEventID() << endl;
    
    G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
    G4int n_trajectories = 0;
    if(trajectoryContainer){ 
      n_trajectories = trajectoryContainer->entries(); 
    }
    //  G4cout << "    " << n_trajectories 
    //	 << " trajectories stored in this event." << endl;
    
    if(G4VVisManager::GetConcreteInstance()){
      for(G4int i=0; i<n_trajectories; i++) {
	G4Trajectory* trj = (G4Trajectory*)(*(evt->GetTrajectoryContainer()))[i];
	if (drawFlag == "all") trj->DrawTrajectory(50);
	else if ((drawFlag == "charged")&&(trj->GetCharge() != 0.))
	  trj->DrawTrajectory(50); 
	trj->ShowTrajectory(); 
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

































