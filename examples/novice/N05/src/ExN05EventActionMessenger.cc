// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN05EventActionMessenger.cc,v 1.2 1999-12-15 14:49:30 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "ExN05EventActionMessenger.hh"
#include "ExN05EventAction.hh"

#include "globals.hh"
#include "G4UIcmdWithABool.hh"

ExN05EventActionMessenger::ExN05EventActionMessenger(ExN05EventAction* EA)
:EventAction(EA)
{
  // uses "/event" directory already
  // defined in G4EvManMessenger
  drawEventCmd = new G4UIcmdWithABool("/event/Draw",this);
  drawEventCmd->SetGuidance("Set drawFlag to Draw an event. (default = false)");
  drawEventCmd->SetParameterName("Draw", false);
  drawEventCmd->SetDefaultValue(false);
}

void ExN05EventActionMessenger::SetNewValue(G4UIcommand* command, G4String newValues)
{
  if( command->GetCommandName() == "Draw" )
  {
    G4int vl;
    const char* t = newValues;
    G4std::istrstream is((char*)t);
    is >> vl;
    EventAction->SetDrawFlag(vl!=0);
  }
}

