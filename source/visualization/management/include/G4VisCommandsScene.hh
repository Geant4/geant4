// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsScene.hh,v 1.4 1999-04-15 15:26:26 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/scene commands - John Allison  9th August 1998

#ifndef G4VISCOMMANDSSCENE_HH
#define G4VISCOMMANDSSCENE_HH

#include "G4VVisCommand.hh"

class G4UIcommand;
class G4UIcmdWithAString;
class G4UIcmdWithoutParameter;

class G4VVisCommandScene: public G4VVisCommand {
public:
  // Uses compiler defaults for destructor, copy constructor and assignment.
  G4VVisCommandScene ();
protected:
  void UpdateCandidateLists ();
  void UpdateVisManagerSceneAndViewParameters (const G4String& sceneName = "");
};

class G4VisCommandSceneCreate: public G4VVisCommandScene {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneCreate ();
  ~G4VisCommandSceneCreate ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4String NextName ();
  G4UIcmdWithAString* fpCommand;
  G4int fId;
};

class G4VisCommandSceneList: public G4VVisCommandScene {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneList ();
  ~G4VisCommandSceneList ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneNotifyHandlers: public G4VVisCommandScene {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneNotifyHandlers ();
  ~G4VisCommandSceneNotifyHandlers ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneRemove: public G4VVisCommandScene {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneRemove ();
  ~G4VisCommandSceneRemove ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneSelect: public G4VVisCommandScene {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneSelect ();
  ~G4VisCommandSceneSelect ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithAString* fpCommand;
};

#endif
