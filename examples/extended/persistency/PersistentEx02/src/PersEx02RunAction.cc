// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: PersEx02RunAction.cc,v 1.3 1999/11/29 18:33:28 morita Exp $
// GEANT4 tag $Name: geant4-01-01 $
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

