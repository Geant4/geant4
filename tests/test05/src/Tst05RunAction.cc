
#include "Tst05RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

Tst05RunAction::Tst05RunAction()
{
  timer = new G4Timer;
}

Tst05RunAction::~Tst05RunAction()
{
  delete timer;
}

void Tst05RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/event/Verbose 1");
  UI->ApplyCommand("/tracking/Verbose 1");

  //G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  timer->Start();
}

void Tst05RunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();
  //G4cout << "number of event = " << aRun->GetNumberOfEvent() 
  //     << " " << *timer << G4endl;
}

