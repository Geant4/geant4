
#include "Tst01RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

Tst01RunAction::Tst01RunAction()
{
  runIDcounter = 0;
}

Tst01RunAction::~Tst01RunAction()
{
}

void Tst01RunAction::BeginOfRunAction(G4Run* aRun)
{
  aRun->SetRunID(runIDcounter++);
}

void Tst01RunAction::EndOfRunAction(G4Run* aRun)
{
}

