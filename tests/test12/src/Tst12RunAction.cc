// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst12RunAction.cc,v 1.2 1999-04-17 08:34:18 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst12RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

Tst12RunAction::Tst12RunAction()
{
  runIDcounter = 0;
}

Tst12RunAction::~Tst12RunAction()
{
}

void Tst12RunAction::BeginOfRunAction(const G4Run* aRun)
{
  ((G4Run*)(aRun))->SetRunID(runIDcounter++);
}

void Tst12RunAction::EndOfRunAction(const G4Run*)
{
}

