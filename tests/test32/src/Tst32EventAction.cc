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
// $Id: Tst32EventAction.cc,v 1.2 2006-06-29 21:59:13 gunter Exp $
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



