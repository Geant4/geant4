// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05RunAction.cc,v 1.3 1999-04-22 21:45:27 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

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
  G4cout << "### Run " << aRun->GetRunID() << " start." << endl;
  //  timer->Start();
}

void ExN05RunAction::EndOfRunAction(const G4Run* aRun)
{
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() << endl;
  //       << " " << *timer << endl;
}

