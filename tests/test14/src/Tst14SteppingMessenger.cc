// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst14SteppingMessenger.cc,v 1.1 1999-05-29 14:12:14 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
#include "Tst14SteppingMessenger.hh"

#include "Tst14SteppingAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "globals.hh"

Tst14SteppingMessenger::Tst14SteppingMessenger(Tst14SteppingAction* SA)
:steppingAction (SA)
{
  steppingDir = new G4UIdirectory("/stepping/");
  steppingDir->SetGuidance("stepping control");
}

Tst14SteppingMessenger::~Tst14SteppingMessenger()
{
  delete steppingDir;
}

void Tst14SteppingMessenger::
         SetNewValue(G4UIcommand* command,G4String newValues)
{

}

   
