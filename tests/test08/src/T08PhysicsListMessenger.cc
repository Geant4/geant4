// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08PhysicsListMessenger.cc,v 1.1 1999-01-08 16:35:20 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "T08PhysicsListMessenger.hh"

#include "T08PhysicsList.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "globals.hh"

T08PhysicsListMessenger::T08PhysicsListMessenger(T08PhysicsList * List)
:myList(List)
{
  ProcessCmd = new G4UIcmdWithoutParameter("/run/setEmProcess",this);
  ProcessCmd->SetGuidance("select electromagnetic processes");
}

T08PhysicsListMessenger::~T08PhysicsListMessenger()
{
  delete ProcessCmd;
}
  
void T08PhysicsListMessenger::SetNewValue(G4UIcommand* command,G4String newValues)
{
  if(command == ProcessCmd) myList->SetStatusEmProcess();
}

