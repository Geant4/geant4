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
// $Id: G4VisCommandsScene.hh,v 1.17 2006-06-29 21:28:40 gunter Exp $
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
  G4VVisCommandScene ();
  virtual ~G4VVisCommandScene ();
protected:
  G4String CurrentSceneName ();
private:
  G4VVisCommandScene (const G4VVisCommandScene&);
  G4VVisCommandScene& operator = (const G4VVisCommandScene&);
};

class G4VisCommandSceneCreate: public G4VVisCommandScene {
public:
  G4VisCommandSceneCreate ();
  virtual ~G4VisCommandSceneCreate ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneCreate (const G4VisCommandSceneCreate&);
  G4VisCommandSceneCreate& operator = (const G4VisCommandSceneCreate&);
  G4String NextName ();
  G4UIcmdWithAString* fpCommand;
  G4int fId;
};

class G4VisCommandSceneEndOfEventAction: public G4VVisCommandScene {
public:
  G4VisCommandSceneEndOfEventAction ();
  virtual ~G4VisCommandSceneEndOfEventAction ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneEndOfEventAction (const G4VisCommandSceneEndOfEventAction&);
  G4VisCommandSceneEndOfEventAction& operator =
  (const G4VisCommandSceneEndOfEventAction&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneEndOfRunAction: public G4VVisCommandScene {
public:
  G4VisCommandSceneEndOfRunAction ();
  virtual ~G4VisCommandSceneEndOfRunAction ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneEndOfRunAction (const G4VisCommandSceneEndOfRunAction&);
  G4VisCommandSceneEndOfRunAction& operator =
  (const G4VisCommandSceneEndOfRunAction&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneList: public G4VVisCommandScene {
public:
  G4VisCommandSceneList ();
  virtual ~G4VisCommandSceneList ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneList (const G4VisCommandSceneList&);
  G4VisCommandSceneList& operator = (const G4VisCommandSceneList&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneNotifyHandlers: public G4VVisCommandScene {
public:
  G4VisCommandSceneNotifyHandlers ();
  virtual ~G4VisCommandSceneNotifyHandlers ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneNotifyHandlers (const G4VisCommandSceneNotifyHandlers&);
  G4VisCommandSceneNotifyHandlers& operator =
  (const G4VisCommandSceneNotifyHandlers&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneSelect: public G4VVisCommandScene {
public:
  G4VisCommandSceneSelect ();
  virtual ~G4VisCommandSceneSelect ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneSelect (const G4VisCommandSceneSelect&);
  G4VisCommandSceneSelect& operator = (const G4VisCommandSceneSelect&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneTransientsAction: public G4VVisCommandScene {
public:
  G4VisCommandSceneTransientsAction ();
  virtual ~G4VisCommandSceneTransientsAction ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneTransientsAction (const G4VisCommandSceneTransientsAction&);
  G4VisCommandSceneTransientsAction& operator =
  (const G4VisCommandSceneTransientsAction&);
  G4UIcmdWithAString* fpCommand;
};

#endif
