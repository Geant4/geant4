// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst02RunAction.cc,v 1.2 1999-04-16 09:41:28 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst02RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

Tst02RunAction::Tst02RunAction()
{
  runIDcounter = 0;
}

Tst02RunAction::~Tst02RunAction()
{
}

void Tst02RunAction::BeginOfRunAction(const G4Run* aRun)
{
  ((G4Run*)(aRun))->SetRunID(runIDcounter++);
}

void Tst02RunAction::EndOfRunAction(const G4Run* aRun)
{
}

