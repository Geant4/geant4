// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsSceneHandler.hh,v 1.2 1999-01-09 16:31:02 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/sceneHandler commands - John Allison  10th October 1998

#ifndef G4VISCOMMANDSSCENEHANDLER_HH
#define G4VISCOMMANDSSCENEHANDLER_HH

#include "G4VVisCommand.hh"

class G4UIcommand;
class G4UIcmdWithAString;

class G4VVisCommandSceneHandler: public G4VVisCommand {
public:
  // Uses compiler defaults for destructor, copy constructor and assignment.
  G4VVisCommandSceneHandler ();
protected:
  void UpdateCandidateLists ();
  static G4String fSceneHandlerNameList;
  // member so that it has long life - static because shared between objects.
};

class G4VisCommandSceneHandlerAttach: public G4VVisCommandSceneHandler {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneHandlerAttach ();
  ~G4VisCommandSceneHandlerAttach ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneHandlerCreate: public G4VVisCommandSceneHandler {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneHandlerCreate ();
  ~G4VisCommandSceneHandlerCreate ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4String NextName ();
  G4UIcommand* fpCommand;
  G4int fId;
};

class G4VisCommandSceneHandlerList: public G4VVisCommandSceneHandler {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneHandlerList ();
  ~G4VisCommandSceneHandlerList ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneHandlerRemove: public G4VVisCommandSceneHandler {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneHandlerRemove ();
  ~G4VisCommandSceneHandlerRemove ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneHandlerSelect: public G4VVisCommandSceneHandler {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneHandlerSelect ();
  ~G4VisCommandSceneHandlerSelect ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithAString* fpCommand;
};

#endif
