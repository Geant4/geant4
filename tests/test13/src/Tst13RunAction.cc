// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst13RunAction.cc,v 1.1 1999-01-08 16:35:47 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst13RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

Tst13RunAction::Tst13RunAction()
{
  runIDcounter = 0;
}

Tst13RunAction::~Tst13RunAction()
{
}

void Tst13RunAction::BeginOfRunAction(G4Run* aRun)
{
  aRun->SetRunID(runIDcounter++);
}

void Tst13RunAction::EndOfRunAction(G4Run* aRun)
{
}

