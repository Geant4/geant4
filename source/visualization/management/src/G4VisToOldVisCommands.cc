// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisToOldVisCommands.cc,v 1.2 1999-01-09 16:31:35 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Implements some /vis/ commands as /vis~/ temporarily.
// John Allison  9th December 1998.

#include "G4VisToOldVisCommands.hh"

#include "G4UImanager.hh"
#include "G4UIcommandTree.hh"
#include "G4UIdirectory.hh"

G4VisToOldVisCommands::G4VisToOldVisCommands () {

  G4UIcommandTree* commandTree = G4UImanager::GetUIpointer () -> GetTree ();
  G4UIcommandTree* oldVisCommandTree = commandTree -> GetTree ("/vis~/");
  G4int nOldVisTreeEntries = oldVisCommandTree -> GetTreeEntry ();
  int i;
  for (i = 1; i <= nOldVisTreeEntries; i++) {
    G4UIcommandTree* subTree = oldVisCommandTree -> GetTree (i);
    G4String subTreeName = subTree -> GetPathName ();
    G4String newDirectoryName = subTreeName.replace (0, 6, "/vis/");
    newDirectoryName = newDirectoryName (0, newDirectoryName.last ('/') + 1);
    // G4cout << "NewDirectoryName: " << newDirectoryName << endl;
    new G4UIdirectory (newDirectoryName);
    G4int nCommands = subTree -> GetCommandEntry ();
    int j;
    for (j = 1; j <= nCommands; j++) {
      G4UIcommand* command = subTree -> GetCommand (j);
      G4String commandPath = command -> GetCommandPath ();
      // G4cout << "SubCommandPath:   " << commandPath << endl;
      G4String newCommandPath = commandPath.replace (0, 6, "/vis/");
      // G4cout << "NewSubCommandPath: " << newCommandPath << endl;
      G4UIcommand* newCommand = new G4UIcommand (newCommandPath, this);
      newCommand -> SetRange (command -> GetRange ());
      G4int nGuidances = command -> GetGuidanceEntries ();
      int k;
      for (k = 0; k < nGuidances; k++) {
	newCommand -> SetGuidance (command -> GetGuidanceLine (k));
      }
      G4int nParameters = command -> GetParameterEntries ();
      for (k = 0; k < nParameters; k++) {
	newCommand -> SetParameter (command -> GetParameter (k));
      }
    }
  }
}

G4VisToOldVisCommands::~G4VisToOldVisCommands () {}

void G4VisToOldVisCommands::SetNewValue
(G4UIcommand* command, G4String newValues) {
  G4String commandPath = command -> GetCommandPath ();
  if (commandPath == "/vis/create_view/new_graphics_system") {
    G4cout <<
      "Command superceded.  Please use:"
      "\n  /vis/scene/create"
      "\n  /vis/sceneHandler/create " << newValues <<
      "\n  /vis/viewer/create"
      "\nThis way you can also use \"/vis/scene/addVolume\" to select"
      "\n  sub-detector components.  (The command is still available"
      "\n  as \"/vis~/create_view/new_graphics/system " << newValues <<
      "\".)"
	   << endl;
  }
  else {
    G4String oldCommandPath = commandPath.replace (0, 5, "/vis~/");
    G4UImanager::GetUIpointer () ->
      ApplyCommand (oldCommandPath + " " + newValues);
  }
}
