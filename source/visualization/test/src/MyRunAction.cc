// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyRunAction.cc,v 1.1 1999-04-16 10:32:37 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MyRunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

#ifdef G4VIS_USE_OPACS
OHistogram MyRunAction::myHisto1D = NULL;
OHistogram MyRunAction::myHisto2D = NULL;
#endif

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
  //aRun->transient(true);

  G4UImanager* UI = G4UImanager::GetUIpointer();
  UI->ApplyCommand("/tracking/storeTrajectory 1");

#ifdef G4VIS_USE_OPACS
  // Histo are created in an osh script somewhere (test19.odb).
  myHisto1D = OHistogramGetIdentifier("HepJamesRandom1");
  myHisto2D = OHistogramGetIdentifier("HepJamesVsDRand48");
#endif

  G4cout << "### Run " << aRun->GetRunID() << " start." << endl;
  timer->Start();
}

void MyRunAction::EndOfRunAction(G4Run* aRun)
{
  timer->Stop();
  G4cout << "number of event = " << aRun->GetNumberOfEvent() 
       << " " << *timer << endl;
#ifdef G4VIS_USE_OPACS
  myHisto1D = NULL;
  myHisto2D = NULL;
#endif
}

