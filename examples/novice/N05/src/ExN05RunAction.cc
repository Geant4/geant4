// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05RunAction.cc,v 1.5 2000-01-06 15:06:51 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// Make this appear first!
#include "G4Timer.hh"

#include "ExN05RunAction.hh"
#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

ExN05RunAction::ExN05RunAction()
{
  timer = new G4Timer;
}

ExN05RunAction::~ExN05RunAction()
{
  delete timer;
}

void ExN05RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  //  timer->Start();
}

void ExN05RunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() << G4endl;
  //       << " " << *timer << G4endl;
}

