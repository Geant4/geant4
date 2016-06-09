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
// $Id: G4VisCommandsSceneAdd.hh,v 1.21 2010-06-03 10:17:44 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// /vis/scene commands - John Allison  9th August 1998

#ifndef G4VISCOMMANDSSCENEADD_HH
#define G4VISCOMMANDSSCENEADD_HH

#include "G4VisCommandsScene.hh"

class G4UIcmdWithoutParameter;
class G4UIcmdWithAString;

#include "G4Transform3D.hh"
#include "G4VisAttributes.hh"

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

class G4VisCommandSceneAddDigis: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddDigis ();
  virtual ~G4VisCommandSceneAddDigis ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddDigis (const G4VisCommandSceneAddDigis&);
  G4VisCommandSceneAddDigis& operator = (const G4VisCommandSceneAddDigis&);
  G4UIcmdWithoutParameter* fpCommand;
};

class G4VisCommandSceneAddEventID: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddEventID ();
  virtual ~G4VisCommandSceneAddEventID ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddEventID (const G4VisCommandSceneAddEventID&);
  G4VisCommandSceneAddEventID& operator = (const G4VisCommandSceneAddEventID&);
  struct EventID {
    EventID(G4VisManager* vm, G4int size, G4double x, G4double y):
      fpVisManager(vm), fSize(size), fX(x), fY(y) {}
    void operator()(G4VGraphicsScene&, const G4Transform3D&);
    G4VisManager* fpVisManager;
    G4int fSize;
    G4double fX, fY;
  };
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

class G4VisCommandSceneAddLogo: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddLogo ();
  virtual ~G4VisCommandSceneAddLogo ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddLogo (const G4VisCommandSceneAddLogo&);
  G4VisCommandSceneAddLogo& operator = (const G4VisCommandSceneAddLogo&);
  class G4Logo {
  public:
    G4Logo(G4double height, const G4VisAttributes&);
    ~G4Logo();
    void operator()(G4VGraphicsScene&, const G4Transform3D&);
  private:
    G4double fHeight;
    G4VisAttributes fVisAtts;
    G4Polyhedron *fpG, *fp4;
  };
  G4UIcommand* fpCommand;
};

class G4VisCommandSceneAddPSHits: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddPSHits ();
  virtual ~G4VisCommandSceneAddPSHits ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddPSHits (const G4VisCommandSceneAddPSHits&);
  G4VisCommandSceneAddPSHits& operator = (const G4VisCommandSceneAddPSHits&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneAddScale: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddScale ();
  virtual ~G4VisCommandSceneAddScale ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddScale (const G4VisCommandSceneAddScale&);
  G4VisCommandSceneAddScale& operator = (const G4VisCommandSceneAddScale&);
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
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSceneAddUserAction: public G4VVisCommandScene {
public:
  G4VisCommandSceneAddUserAction ();
  virtual ~G4VisCommandSceneAddUserAction ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSceneAddUserAction (const G4VisCommandSceneAddUserAction&);
  G4VisCommandSceneAddUserAction& operator = (const G4VisCommandSceneAddUserAction&);
  G4UIcommand* fpCommand;
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
