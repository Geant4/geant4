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
// $Id: Mars01EventAction.cc,v 1.1 2001-12-13 14:58:47 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Mars01EventAction.hh"
#include "Mars01Hit.hh"

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

Mars01EventAction::Mars01EventAction()
{
  TallyFile.open("mars.out");
}

Mars01EventAction::~Mars01EventAction()
{
  TallyFile.close();
}

void Mars01EventAction::BeginOfEventAction(const G4Event*)
{

}

void Mars01EventAction::EndOfEventAction(const G4Event* evt )
{
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  G4int collID    = SDman->GetCollectionID(colNam="HitsCollection");
  
  G4int nevt = evt->GetEventID();
  
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  Mars01HitsCollection* pHC    = 0;
  if(HCE) {
    pHC    = (Mars01HitsCollection*)(HCE->GetHC(collID));
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



