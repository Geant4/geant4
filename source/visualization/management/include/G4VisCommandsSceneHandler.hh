//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4VisCommandsSceneHandler.hh,v 1.6 2001-07-11 10:09:16 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/sceneHandler commands - John Allison  10th October 1998

#ifndef G4VISCOMMANDSSCENEHANDLER_HH
#define G4VISCOMMANDSSCENEHANDLER_HH

#include "G4VVisCommand.hh"

class G4UIcommand;
class G4UIcmdWithAString;

class G4VVisCommandSceneHandler: public G4VVisCommand {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VVisCommandSceneHandler ();
  ~G4VVisCommandSceneHandler ();
protected:
  void UpdateCandidateLists ();
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
