// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst17SteppingMessenger.cc,v 1.2 1999-12-15 14:54:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
#include "Tst17SteppingMessenger.hh"

#include "Tst17SteppingAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "globals.hh"

Tst17SteppingMessenger::Tst17SteppingMessenger(Tst17SteppingAction* SA)
:steppingAction (SA)
{
  steppingDir = new G4UIdirectory("/stepping/");
  steppingDir->SetGuidance("stepping control");
}

Tst17SteppingMessenger::~Tst17SteppingMessenger()
{
  delete steppingDir;
}

void Tst17SteppingMessenger::
         SetNewValue(G4UIcommand* command,G4String newValues)
{

}

   
