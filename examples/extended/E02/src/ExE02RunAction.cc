
#include "ExE02RunAction.hh"

#include "G4Run.hh"
#include "G4VPersistencyManager.hh"
#include "G4TransportationManager.hh"
#include "G4VPhysicalVolume.hh"

ExE02RunAction::ExE02RunAction()
{
  runIDcounter = 0;
}

ExE02RunAction::~ExE02RunAction()
{
}

void ExE02RunAction::BeginOfRunAction(G4Run* aRun)
{
  aRun->SetRunID(runIDcounter++);
  if(runIDcounter==1)
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

void ExE02RunAction::EndOfRunAction(G4Run* aRun)
{
}

