// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

#include "Tst03RunAction.hh"

#include "G4ios.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"

Tst03RunAction::Tst03RunAction()
{
  timer = new G4Timer;
  runIDcounter = 0;
}

Tst03RunAction::~Tst03RunAction()
{
  delete timer;
}

void Tst03RunAction::BeginOfRunAction(G4Run* aRun)
{
  aRun->SetRunID(runIDcounter++);

  G4UImanager* UI = G4UImanager::GetUIpointer();
  //UI->ApplyCommand("/event/verbose 1");
  //UI->ApplyCommand("/tracking/verbose 1");

  G4cout << "### Run " << aRun->GetRunID() << " start." << endl;
  timer->Start();
}

void Tst03RunAction::EndOfRunAction(G4Run* aRun)
{
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() 
       << " " << *timer << endl;
}
