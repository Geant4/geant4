#include "ExGflashRunAction.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"




ExGflashRunAction::ExGflashRunAction(): runID(0)
{
}


ExGflashRunAction::~ExGflashRunAction()
{
}


void ExGflashRunAction::BeginOfRunAction(const G4Run* aRun)
{  
   ((G4Run *)(aRun))->SetRunID(runID++);
   
   std::cout << "### Run " << aRun->GetRunID() << " start." << std::endl;
   
   if (G4VVisManager::GetConcreteInstance())
     {
       G4UImanager* UI = G4UImanager::GetUIpointer();
       UI->ApplyCommand("/vis/scene/notifyHandlers");
     } 
}


void ExGflashRunAction::EndOfRunAction(const G4Run* aRun)
{ 
  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
    }

  G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;

}











