// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MySteppingActionMessenger.cc,v 1.1 1999-01-07 16:05:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "MySteppingActionMessenger.hh"

#include "MySteppingAction.hh"
#include "G4UIparameter.hh"
#include "globals.hh"

MySteppingActionMessenger::MySteppingActionMessenger(MySteppingAction * mySA)
:mySteppingAction(mySA)
{
  G4UIcommand * command;
  G4UIparameter * param;

  command = new G4UIcommand("/Step/Draw",this);
  command->SetGuidance("Set drawFlag to Draw event Step.");
  command->SetGuidance("  default value is 0 (false).");
  param = new G4UIparameter("flag (1:true/0:fasle)",'i',true);
  param->SetDefaultValue(0);
  command->SetParameter(param);
  AddUIcommand(command);
}

void MySteppingActionMessenger::SetNewValue(G4UIcommand * command,G4String newValues)
{
  if( command->GetCommandName() == "Draw" )
  {
    G4int vl;
    const char* t = newValues;
    istrstream is((char*)t);
    is >> vl;
    mySteppingAction->SetDrawFlag(vl!=0);
  }
}

