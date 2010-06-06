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
//    **************************************
//    *                                    *
//    *        CellEventAction.cc          *
//    *                                    *
//    **************************************
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//	   Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
//
// History:
// -----------
// 20 September 2006   S. Guatelli, B. Mascialino   1st implementation
//
// -------------------------------------------------------------------
 
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ios.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#ifdef G4ANALYSIS_USE
#include "CellAnalysisManager.hh"
#endif
#include "CellTrackerHit.hh"
#include "CellEventAction.hh"

CellEventAction::CellEventAction()
{ }

CellEventAction::~CellEventAction()
{ }

void CellEventAction::BeginOfEventAction(const G4Event*)
{ 
if (collisionID==-1)
    {
      G4SDManager * SDman = G4SDManager::GetSDMpointer();
      collisionID = SDman->GetCollectionID("TstCellCollection");
    }  

 totalEnergy = 0;
}
 
void CellEventAction::EndOfEventAction(const G4Event* evt)
{
}

G4int CellEventAction::GetEventNo()
{
  G4int evno = fpEventManager -> GetConstCurrentEvent() -> GetEventID() ;
  return evno ;
}

