// $Id: Tst10RunAction.cc,v 1.1 1999-01-08 16:35:34 gunter Exp $
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

void Tst10RunAction::BeginOfRunAction(G4Run* aRun)
{
  aRun->SetRunID(runIDcounter++);
}

void Tst10RunAction::EndOfRunAction(G4Run* aRun)
{
}

