#include "Tst34RunAction.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"




Tst34RunAction::Tst34RunAction(): runID(0)
{
}


Tst34RunAction::~Tst34RunAction()
{
}


void Tst34RunAction::BeginOfRunAction(const G4Run* aRun)
{  
   ((G4Run *)(aRun))->SetRunID(runID++);
   
   std::cout << "### Run " << aRun->GetRunID() << " start." << std::endl;
}


void Tst34RunAction::EndOfRunAction(const G4Run* aRun)
{ 

  G4cout << "number of events = " << aRun->GetNumberOfEvent() << G4endl;

}











