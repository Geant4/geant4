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
// $Id: G4VisCommandsScene.hh,v 1.12 2001-09-10 10:49:18 johna Exp $
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
  void UpdateCandidateLists ();
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

class G4VisCommandSceneEdit: public G4VVisCommandScene {
public:
  G4VisCommandSceneEdit ();
  virtual ~G4VisCommandSceneEdit ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneEdit (const G4VisCommandSceneEdit&);
  G4VisCommandSceneEdit& operator = (const G4VisCommandSceneEdit&);
  G4UIcmdWithAString* fpCommand;
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
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneRemove: public G4VVisCommandScene {
public:
  G4VisCommandSceneRemove ();
  virtual ~G4VisCommandSceneRemove ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneRemove (const G4VisCommandSceneRemove&);
  G4VisCommandSceneRemove& operator = (const G4VisCommandSceneRemove&);
  G4UIcmdWithAString* fpCommand;
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

#endif
