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
// $Id: PersEx02RunAction.cc,v 1.3.8.1 2001/06/28 19:07:25 gunter Exp $
// GEANT4 tag $Name:  $
//

#include "PersEx02RunAction.hh"

#include "G4Run.hh"
#include "G4PersistencyManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"

PersEx02RunAction::PersEx02RunAction()
{
}

PersEx02RunAction::~PersEx02RunAction()
{
}

void PersEx02RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4PersistencyManager* persM 
    = G4PersistencyManager::GetPersistencyManager();

  if(aRun->GetRunID()==0)
  {
    if(persM)
    {
      G4VPhysicalVolume* theWorld
        = G4TransportationManager::GetTransportationManager()
          ->GetNavigatorForTracking()->GetWorldVolume();
      persM->Store(theWorld);
    }
  }

  // start event DB transaction with sustained mode
  persM->StartTransaction(kEventDB, kUpdate, true);

}

void PersEx02RunAction::EndOfRunAction(const G4Run* )
{
  G4PersistencyManager* persM 
    = G4PersistencyManager::GetPersistencyManager();
  if(persM)
  {
    // commit the sustained event DB transaction
    persM->Commit(kEventDB, true);
  }
}

