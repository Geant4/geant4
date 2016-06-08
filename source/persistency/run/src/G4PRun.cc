// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PRun.cc,v 1.5 1999/12/15 14:51:27 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//

#include "G4Run.hh"
#include "G4PRun.hh"
#include "G4Event.hh"

G4PRun::G4PRun()
// :runID(0),numberOfEvent(0),HCtable(NULL),DCtable(NULL)
:runID(0),numberOfEvent(0)
{;}

G4PRun::G4PRun(const G4Run* aRun)
{
  runID = aRun->GetRunID();
  numberOfEvent = aRun->GetNumberOfEvent();
}

G4PRun::~G4PRun()
{;}

G4Run* G4PRun::MakeTransientObject()
{
  G4Run* aRun = new G4Run();

  aRun->SetRunID(runID);

  G4Event* dummyEvt = new G4Event;

  for(G4int i = 0; i < numberOfEvent; i++)
  {
    aRun->RecordEvent(dummyEvt);
  }

  delete dummyEvt;

  return aRun;
}

