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
// $Id: G4VisCommandsSceneAdd.hh,v 1.8 2001-07-11 10:09:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/scene commands - John Allison  9th August 1998

#ifndef G4VISCOMMANDSSCENEADD_HH
#define G4VISCOMMANDSSCENEADD_HH

#include "G4VisCommandsScene.hh"

class G4VisCommandSceneAddAxes: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddAxes ();
  virtual ~G4VisCommandSceneAddAxes ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddAxes (const G4VisCommandSceneAddAxes&);
  G4VisCommandSceneAddAxes& operator = (const G4VisCommandSceneAddAxes&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddGhosts: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddGhosts ();
  virtual ~G4VisCommandSceneAddGhosts ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddGhosts (const G4VisCommandSceneAddGhosts&);
  G4VisCommandSceneAddGhosts& operator =
  (const G4VisCommandSceneAddGhosts&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneAddHits: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddHits ();
  virtual ~G4VisCommandSceneAddHits ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddHits (const G4VisCommandSceneAddHits&);
  G4VisCommandSceneAddHits& operator = (const G4VisCommandSceneAddHits&);
  G4UIcmdWithoutParameter* fpCommand;
};

class G4VisCommandSceneAddLogicalVolume: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddLogicalVolume ();
  virtual ~G4VisCommandSceneAddLogicalVolume ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddLogicalVolume (const G4VisCommandSceneAddLogicalVolume&);
  G4VisCommandSceneAddLogicalVolume& operator =
  (const G4VisCommandSceneAddLogicalVolume&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddText: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddText ();
  virtual ~G4VisCommandSceneAddText ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddText (const G4VisCommandSceneAddText&);
  G4VisCommandSceneAddText& operator = (const G4VisCommandSceneAddText&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddTrajectories: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddTrajectories ();
  virtual ~G4VisCommandSceneAddTrajectories ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddTrajectories (const G4VisCommandSceneAddTrajectories&);
  G4VisCommandSceneAddTrajectories& operator =
  (const G4VisCommandSceneAddTrajectories&);
  G4UIcmdWithoutParameter* fpCommand;
};

class G4VisCommandSceneAddVolume: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddVolume ();
  virtual ~G4VisCommandSceneAddVolume ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddVolume (const G4VisCommandSceneAddVolume&);
  G4VisCommandSceneAddVolume& operator = (const G4VisCommandSceneAddVolume&);
  G4UIcommand* fpCommand;
};

#endif
