// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsSceneInclude.hh,v 1.2 1999-01-09 16:31:03 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/scene commands - John Allison  9th August 1998

#ifndef G4VISCOMMANDSSCENEINCLUDE_HH
#define G4VISCOMMANDSSCENEINCLUDE_HH

#include "G4VisCommandsScene.hh"

class G4VisCommandSceneIncludeHits: public G4VVisCommandScene {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneIncludeHits ();
  ~G4VisCommandSceneIncludeHits ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithoutParameter* fpCommand;
};

class G4VisCommandSceneIncludeTrajectories: public G4VVisCommandScene {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneIncludeTrajectories ();
  ~G4VisCommandSceneIncludeTrajectories ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithoutParameter* fpCommand;
};

#endif
