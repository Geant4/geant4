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
// $Id: Tst32EventAction.cc,v 1.1 2002-06-13 12:16:35 jwellisc Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst32EventAction.hh"
#include "Tst32Hit.hh"

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

Tst32EventAction::Tst32EventAction()
{
  TallyFile.open("mars.out");
}

Tst32EventAction::~Tst32EventAction()
{
  TallyFile.close();
}

void Tst32EventAction::BeginOfEventAction(const G4Event*)
{

}

void Tst32EventAction::EndOfEventAction(const G4Event* evt )
{
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  G4int collID    = SDman->GetCollectionID(colNam="HitsCollection");
  
  G4int nevt = evt->GetEventID();
  
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  Tst32HitsCollection* pHC    = 0;
  if(HCE) {
    pHC    = (Tst32HitsCollection*)(HCE->GetHC(collID));
  }
  
  if(pHC){
    int n_hit = pHC->entries();
    TallyFile << nevt << ' ';
    for(int i=0;i<n_hit;i++){ 
      G4double totE = (*pHC)[i]->GetEdep(); 
      TallyFile << totE/MeV << ' ';  
    }
    TallyFile << G4endl;
  }
}



