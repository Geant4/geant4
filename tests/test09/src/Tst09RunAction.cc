// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst09RunAction.cc,v 1.2 1999-04-17 07:47:03 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst09RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

Tst09RunAction::Tst09RunAction()
{
  runIDcounter = 0;
}

Tst09RunAction::~Tst09RunAction()
{
}

void Tst09RunAction::BeginOfRunAction(const G4Run* aRun)
{
  ((G4Run*)(aRun))->SetRunID(runIDcounter++);
}

void Tst09RunAction::EndOfRunAction(const G4Run* )
{
}

