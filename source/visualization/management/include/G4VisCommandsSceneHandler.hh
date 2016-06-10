//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4VisCommandsSceneHandler.hh 66373 2012-12-18 09:41:34Z gcosmo $

// /vis/sceneHandler commands - John Allison  10th October 1998

#ifndef G4VISCOMMANDSSCENEHANDLER_HH
#define G4VISCOMMANDSSCENEHANDLER_HH

#include "G4VVisCommand.hh"

class G4UIcommand;
class G4UIcmdWithAString;

class G4VisCommandSceneHandlerAttach: public G4VVisCommand {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneHandlerAttach ();
  ~G4VisCommandSceneHandlerAttach ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneHandlerCreate: public G4VVisCommand {
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

class G4VisCommandSceneHandlerList: public G4VVisCommand {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSceneHandlerList ();
  ~G4VisCommandSceneHandlerList ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneHandlerSelect: public G4VVisCommand {
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
