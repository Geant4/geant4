
#include "MyRunAction.hh"

#include "G4Run.hh"
#include "G4VPersistencyManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"

MyRunAction::MyRunAction()
{
}

MyRunAction::~MyRunAction()
{
}

void MyRunAction::BeginOfRunAction(const G4Run* aRun)
{
  if(aRun->GetRunID()==0)
  {
    G4VPersistencyManager* persM 
      = G4VPersistencyManager::GetPersistencyManager();
    if(persM)
    {
      G4VPhysicalVolume* theWorld
        = G4TransportationManager::GetTransportationManager()
          ->GetNavigatorForTracking()->GetWorldVolume();
      persM->Store(theWorld);
    }
  }
}

void MyRunAction::EndOfRunAction(const G4Run*)
{
}

