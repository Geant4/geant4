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
// $Id: T07EventAction.cc,v 1.9 2006-06-29 21:37:22 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "T07EventAction.hh"

#include "T07CalorHit.hh"
#include "T07EventActionMessenger.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


T07EventAction::T07EventAction()
 : calorimeterCollID(-1),drawFlag("all"),eventMessenger(0)
{
  eventMessenger = new T07EventActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

T07EventAction::~T07EventAction()
{
  delete eventMessenger;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07EventAction::BeginOfEventAction(const G4Event*)
{  if(calorimeterCollID==-1)
  {
    G4SDManager * SDman = G4SDManager::GetSDMpointer();
    calorimeterCollID = SDman->GetCollectionID("CalCollection");
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07EventAction::EndOfEventAction( const G4Event* evt )
{
 
  // if(evt->GetEventID()%100==0)
  //   G4cout << ">>> Event " << evt->GetEventID() << G4endl;

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  T07CalorHitsCollection* CHC = 0;
  if (HCE)
      CHC = (T07CalorHitsCollection*)(HCE->GetHC(calorimeterCollID));

  if (CHC)
  {
    int n_hit = CHC->entries();
    G4double totEAbs=0, totLAbs=0, totEGap=0, totLGap=0;
    for (int i=0;i<n_hit;i++)
      { totEAbs += (*CHC)[i]->GetEdepAbs(); 
        totLAbs += (*CHC)[i]->GetTrakAbs();
        totEGap += (*CHC)[i]->GetEdepGap(); 
        totLGap += (*CHC)[i]->GetTrakGap();
        
      }
    /*
    G4cout
       << "   Absorber: total energy: " << std::setw(7) << G4BestUnit(totEAbs,"Energy")
       << "       total track length: " << std::setw(7) << G4BestUnit(totLAbs,"Length")
       << G4endl
       << "        Gap: total energy: " << std::setw(7) << G4BestUnit(totEGap,"Energy")
       << "       total track length: " << std::setw(7) << G4BestUnit(totLGap,"Length")
       << G4endl;
    */
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
