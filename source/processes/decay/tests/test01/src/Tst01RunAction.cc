// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst01RunAction.cc,v 1.1 2001-02-08 08:41:51 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "Tst01RunAction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"

Tst01RunAction::Tst01RunAction()
{
}

Tst01RunAction::~Tst01RunAction()
{
}

void Tst01RunAction::BeginOfRunAction(const G4Run* aRun)
{
  ftimer.Start();
}

void Tst01RunAction::EndOfRunAction(const G4Run* aRun)
{
   ftimer.Stop();
   G4cout << ftimer <<G4endl;
}

