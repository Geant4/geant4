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
// $Id: G4VisCommands.hh,v 1.8 2009-03-09 12:42:00 allison Exp $

// /vis/set - John Allison  21st March 2012
// Set quantities for use in appropriate future commands.

#ifndef G4VISCOMMANDSSET_HH
#define G4VISCOMMANDSSET_HH

#include "G4VVisCommand.hh"

class G4UIcommand;
class G4UIcmdWithAnInteger;
class G4UIcmdWithADouble;
class G4UIcmdWithAString;

class G4VisCommandSetArrow3DLineSegmentsPerCircle: public G4VVisCommand {
public:
  G4VisCommandSetArrow3DLineSegmentsPerCircle ();
  virtual ~G4VisCommandSetArrow3DLineSegmentsPerCircle ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSetArrow3DLineSegmentsPerCircle (const G4VisCommandSetArrow3DLineSegmentsPerCircle&);
  G4VisCommandSetArrow3DLineSegmentsPerCircle& operator = (const G4VisCommandSetArrow3DLineSegmentsPerCircle&);
  G4UIcmdWithAnInteger* fpCommand;
};

class G4VisCommandSetColour: public G4VVisCommand {
public:
  G4VisCommandSetColour ();
  virtual ~G4VisCommandSetColour ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSetColour (const G4VisCommandSetColour&);
  G4VisCommandSetColour& operator = (const G4VisCommandSetColour&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSetLineWidth: public G4VVisCommand {
public:
  G4VisCommandSetLineWidth ();
  virtual ~G4VisCommandSetLineWidth ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSetLineWidth (const G4VisCommandSetLineWidth&);
  G4VisCommandSetLineWidth& operator = (const G4VisCommandSetLineWidth&);
  G4UIcmdWithADouble* fpCommand;
};

class G4VisCommandSetTextColour: public G4VVisCommand {
public:
  G4VisCommandSetTextColour ();
  virtual ~G4VisCommandSetTextColour ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSetTextColour (const G4VisCommandSetTextColour&);
  G4VisCommandSetTextColour& operator = (const G4VisCommandSetTextColour&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSetTextLayout: public G4VVisCommand {
public:
  G4VisCommandSetTextLayout ();
  virtual ~G4VisCommandSetTextLayout ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSetTextLayout (const G4VisCommandSetTextLayout&);
  G4VisCommandSetTextLayout& operator = (const G4VisCommandSetTextLayout&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSetTextSize: public G4VVisCommand {
public:
  G4VisCommandSetTextSize ();
  virtual ~G4VisCommandSetTextSize ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSetTextSize (const G4VisCommandSetTextSize&);
  G4VisCommandSetTextSize& operator = (const G4VisCommandSetTextSize&);
  G4UIcmdWithADouble* fpCommand;
};

class G4VisCommandSetTouchable: public G4VVisCommand {
public:
  G4VisCommandSetTouchable ();
  virtual ~G4VisCommandSetTouchable ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSetTouchable (const G4VisCommandSetTouchable&);
  G4VisCommandSetTouchable& operator = (const G4VisCommandSetTouchable&);
  G4UIcommand* fpCommand;
};

#endif
