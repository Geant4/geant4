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
// $Id: Tst23EventAction.cc,v 1.1 2001-12-14 14:53:41 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst23EventAction.hh"
#include "Tst23Hit.hh"

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

Tst23EventAction::Tst23EventAction()
{
}

Tst23EventAction::~Tst23EventAction()
{;}

void Tst23EventAction::BeginOfEventAction(const G4Event*)
{

}

void Tst23EventAction::EndOfEventAction(const G4Event* evt )
{
  G4SDManager * SDman = G4SDManager::GetSDMpointer();
  G4String colNam;
  G4int collID    = SDman->GetCollectionID(colNam="HitsCollection");
  
  G4cout << ">>> Event " << evt->GetEventID() << G4endl;
  
  G4HCofThisEvent * HCE = evt->GetHCofThisEvent();
  Tst23HitsCollection* pHC    = 0;
  if(HCE) {
    pHC    = (Tst23HitsCollection*)(HCE->GetHC(collID));
  }
  
  if(pHC){
    int n_hit = pHC->entries();
    for(int i=0;i<n_hit;i++){ 
      G4double totE = (*pHC)[i]->GetEdep(); 
      G4cout << "     Total energy deposition in "<<i<<"-th BOX: "
	     << totE / GeV << " (GeV)" << G4endl;
    }
  }
}



