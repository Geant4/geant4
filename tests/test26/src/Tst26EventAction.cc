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
//
// $Id: Tst26EventAction.cc,v 1.6 2006-06-29 21:53:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst26EventAction.hh"

#include "Tst26EventMessenger.hh"
#include "Tst26RunAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26EventAction::Tst26EventAction(Tst26RunAction* run)
:Tst26Run(run),
 drawFlag("all"),
 printModulo(100),
 Eth(1.0*keV)
{
  eventMessenger = new Tst26EventMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst26EventAction::~Tst26EventAction()
{
  delete eventMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26EventAction::BeginOfEventAction(const G4Event* evt)
{
  G4int evtNb = evt->GetEventID();     
 
  E1    = 0.0;
  E9    = 0.0;
  E25   = 0.0;
  Eabs1 = 0.0;
  Eabs2 = 0.0;
  Eabs3 = 0.0;
  Eabs4 = 0.0;
  Evert.clear();
  Nvert.clear();

  //printing survey
  if (evtNb%printModulo == 0) 
    G4cout << "\n---> Begin of Event: " << evtNb << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst26EventAction::EndOfEventAction(const G4Event* evt)
{  
  G4int n = Nvert.size();
  G4int nPad = 0;
  if (n > 0) {
    for(G4int i=0; i<n; i++) {
      if (Evert[i] > Eth) nPad++;
    }
  }
  //  G4cout << "EndOfEvent: " << E1 << " " <<  E9 << " " <<  E25 << " " 
  //       <<  Eabs1 << " " << Eabs2 << " " <<  Eabs3 << " " <<  Eabs4 << " " <<  nPad << G4endl;
  Tst26Run->AddEvent(E1, E9, E25, Eabs1, Eabs2, Eabs3, Eabs4, nPad);
  
  // extract the trajectories and draw them
  //
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

void Tst26EventAction::AddEnergy(G4double edep, G4int volIndex, G4int copyNo) 
{  
  if(0 == volIndex) {
    E25 += edep;
    if( (6<=copyNo&&8>=copyNo) || (11<=copyNo&&13>=copyNo) || 
        (16<=copyNo&&18>=copyNo)) {
      E9 += edep;
      if(12 == copyNo) E1 += edep;
    }
  } else if (1 == volIndex) {
    Eabs1 += edep;  
  } else if (2 == volIndex) {
    Eabs2 += edep;  
  } else if (3 == volIndex) {
    Eabs3 += edep;  
  } else if (4 == volIndex) {
    Eabs4 += edep;  
  } else if (5 == volIndex) {
    G4int n = Nvert.size();
    G4bool newPad = true;
    if (n > 0) {
      for(G4int i=0; i<n; i++) {
        if (copyNo == Nvert[i]) {
          newPad = false;
          Evert[i] += edep;
          break;
	}
      }
    }
    if(newPad) {
      Nvert.push_back(copyNo);
      Evert.push_back(edep);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


