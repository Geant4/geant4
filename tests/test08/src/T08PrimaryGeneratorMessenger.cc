// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08PrimaryGeneratorMessenger.cc,v 1.1 1999-01-08 16:35:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "T08PrimaryGeneratorMessenger.hh"

#include "T08PrimaryGeneratorAction.hh"
#include "G4UIcmdWithAString.hh"
#include "globals.hh"

T08PrimaryGeneratorMessenger::T08PrimaryGeneratorMessenger(T08PrimaryGeneratorAction* myGun)
:myAction(myGun)
{ 

  RndmCmd = new G4UIcmdWithAString("/gun/random",this);
  RndmCmd->SetGuidance("Shoot randomly the incident particle.");
  RndmCmd->SetGuidance("  Choice : on, off(default)");
  RndmCmd->SetParameterName("choice",true);
  RndmCmd->SetDefaultValue("off");
  RndmCmd->SetCandidates("on off");
  RndmCmd->AvailableForStates(PreInit,Idle);
}

T08PrimaryGeneratorMessenger::~T08PrimaryGeneratorMessenger()
{
  delete RndmCmd;
}

void T08PrimaryGeneratorMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{ 
  if( command == RndmCmd )
   { myAction->SetRndmFlag(newValues);}
}

