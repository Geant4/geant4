// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyRunAction.cc,v 1.2 1999-12-15 14:48:45 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MyRunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

MyRunAction::MyRunAction()
{
  timer = new G4Timer;
  runIDcounter = 0;
}

MyRunAction::~MyRunAction()
{
  delete timer;
}

void MyRunAction::BeginOfRunAction(G4Run* aRun)
{
  aRun->SetRunID(runIDcounter++);
  aRun->transient(true);

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  timer->Start();
}

void MyRunAction::EndOfRunAction(G4Run* aRun)
{
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() 
       << " " << *timer << G4endl;
}

