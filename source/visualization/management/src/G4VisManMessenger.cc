// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisManMessenger.cc,v 1.1 1999-01-07 16:15:32 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// GEANT4 Visualization Manager Messenger - John Allison 22nd July 1996.

#include "G4VisManMessenger.hh"

#include "G4VisManager.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"

G4VisManMessenger::G4VisManMessenger (G4VisManager* pVMan):
fpVMan (pVMan)
{
  G4UIcommand* command;
  G4UIparameter* param;

  /////////////////////////////////////////////////  /vis~/...  ////
  command = new G4UIcommand ("/vis~/", this);
  command -> SetGuidance ("Deprecated visualization commands.");
  fCommandList.append (command);

  //////////////////////////////////////////////  /vis~/*/...  ////
  AddCommandCamera ();      // Define the /vis~/camera/ sub-commands.
  //  AddCommandClear ();       // Define the /vis~/clear/ sub-commands.
  //  AddCommandCopy ();        // Define the /vis~/copy/ sub-commands.
  //  AddCommandCreateScene (); // Define the /vis~/create_scene/ sub-commands.
  AddCommandCreateView ();  // Define the /vis~/create_view/ sub-commands.
  AddCommandDraw ();        // Define the /vis~/draw/ sub-commands.
  AddCommandLights ();      // Define the /vis~/lights/ sub-commands.
  //  AddCommandPrint ();       // Define the /vis~/print/ sub-commands.
  //  AddCommandRefresh ();     // Define the /vis~/refresh/ sub-commands.
  AddCommandSet ();         // Define the /vis~/set/ sub-commands.
  //  AddCommandShow ();        // Define the /vis~/show/ sub-commands.
  AddCommandExpert ();      // Define the /vis~/expert/ sub-commands.

}

G4VisManMessenger::~G4VisManMessenger () {
  fCommandList.clearAndDestroy ();
}

void G4VisManMessenger::SetNewValue (G4UIcommand* command, G4String newValues)
{
  G4String commandPath = command -> GetCommandPath ();

  if (commandPath.contains ("/vis~/camera/")) {
    DoCommandCamera (commandPath, newValues);
  }

  //  if (commandPath.contains ("/vis~/clear/")) {
  //    DoCommandClear (commandPath, newValues);
  //  }

  //  if (commandPath.contains ("/vis~/copy/")) {
  //    DoCommandCopy (commandPath, newValues);
  //  }

  //  if (commandPath.contains ("/vis~/create_scene/")) {
  //    DoCommandCreateScene (commandPath, newValues);
  //  }

  if (commandPath.contains ("/vis~/create_view/")) {
    DoCommandCreateView (commandPath, newValues);
  }

  if (commandPath.contains ("/vis~/draw/")) {
    DoCommandDraw (commandPath, newValues);
  }

  if (commandPath.contains ("/vis~/lights/")) {
    DoCommandLights (commandPath, newValues);
  }

  //  if (commandPath.contains ("/vis~/print/")) {
  //    DoCommandPrint (commandPath, newValues);
  //  }

  //  if (commandPath.contains ("/vis~/refresh/")) {
  //    DoCommandRefresh (commandPath, newValues);
  //  }

  if (commandPath.contains ("/vis~/set/")) {
    DoCommandSet (commandPath, newValues);
  }

  //  if (commandPath.contains ("/vis~/show/")) {
  //    DoCommandShow (commandPath, newValues);
  //  }

  if (commandPath.contains ("/vis~/expert/")) {
    DoCommandExpert (commandPath, newValues);
  }

}

G4String G4VisManMessenger::GetCurrentValue (G4UIcommand* command) {
  return "";
}

G4bool G4VisManMessenger::ViewValid () {
  if (fpVMan -> IsValidView ()) {
    return true;
  }
  else {
    G4cerr << "Invalid view (have you selected a graphics system?)" << endl;
    return false;
  }
}
