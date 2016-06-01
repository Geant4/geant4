// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05RunAction.cc,v 2.2 1998/07/13 16:34:16 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//

#include "ExN05RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

ExN05RunAction::ExN05RunAction()
{
  timer = new G4Timer;
  runIDcounter = 0;
}

ExN05RunAction::~ExN05RunAction()
{
  delete timer;
}

void ExN05RunAction::BeginOfRunAction(G4Run* aRun)
{
  aRun->SetRunID(runIDcounter++);

  G4cout << "### Run " << aRun->GetRunID() << " start." << endl;
  //  timer->Start();
}

void ExN05RunAction::EndOfRunAction(G4Run* aRun)
{
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() << endl;
  //       << " " << *timer << endl;
}

