
#include "ExE03RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

ExE03RunAction::ExE03RunAction()
{
  timer = new G4Timer;
  runIDcounter = 0;
}

ExE03RunAction::~ExE03RunAction()
{
  delete timer;
}

void ExE03RunAction::BeginOfRunAction(G4Run* aRun)
{
  aRun->SetRunID(runIDcounter++);

  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/event/Verbose 1");
  UI->ApplyCommand("/tracking/Verbose 1");

  G4cout << "### Run " << aRun->GetRunID() << " start." << endl;
  timer->Start();
}

void ExE03RunAction::EndOfRunAction(G4Run* aRun)
{
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() 
       << " " << *timer << endl;
}

