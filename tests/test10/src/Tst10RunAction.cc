// $Id: Tst10RunAction.cc,v 1.2 1999-04-17 08:01:52 kurasige Exp $
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This is a version for maximum particle set
//	History
//        first version              09  Sept. 1998 by S.Magni
// ------------------------------------------------------------

#include "Tst10RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

Tst10RunAction::Tst10RunAction()
{
  runIDcounter = 0;
}

Tst10RunAction::~Tst10RunAction()
{
}

void Tst10RunAction::BeginOfRunAction(const G4Run* aRun)
{
  ((G4Run*)(aRun))->SetRunID(runIDcounter++);
}

void Tst10RunAction::EndOfRunAction(const G4Run* )
{
}

