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

// Compound /vis/ commands - John Allison  15th May 2000

#ifndef G4VISCOMMANDSCOMPOUND_HH
#define G4VISCOMMANDSCOMPOUND_HH

#include "G4VVisCommand.hh"

class G4VisCommandDrawTree: public G4VVisCommand {
public:
  G4VisCommandDrawTree ();
  virtual ~G4VisCommandDrawTree ();
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandDrawTree (const G4VisCommandDrawTree&);
  G4VisCommandDrawTree& operator = (const G4VisCommandDrawTree&);
  G4UIcommand* fpCommand;
};

class G4VisCommandDrawView: public G4VVisCommand {
public:
  G4VisCommandDrawView ();
  virtual ~G4VisCommandDrawView ();
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandDrawView (const G4VisCommandDrawView&);
  G4VisCommandDrawView& operator = (const G4VisCommandDrawView&);
  G4UIcommand* fpCommand;
};

class G4VisCommandDrawLogicalVolume: public G4VVisCommand {
public:
  G4VisCommandDrawLogicalVolume ();
  virtual ~G4VisCommandDrawLogicalVolume ();
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandDrawLogicalVolume (const G4VisCommandDrawLogicalVolume&);
  G4VisCommandDrawLogicalVolume& operator = (const G4VisCommandDrawLogicalVolume&);
  G4UIcommand* fpCommand;
};

  class G4VisCommandDrawVolume: public G4VVisCommand {
public:
  G4VisCommandDrawVolume ();
  virtual ~G4VisCommandDrawVolume ();
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandDrawVolume (const G4VisCommandDrawVolume&);
  G4VisCommandDrawVolume& operator = (const G4VisCommandDrawVolume&);
  G4UIcommand* fpCommand;
};

class G4VisCommandOpen: public G4VVisCommand {
public:
  G4VisCommandOpen ();
  virtual ~G4VisCommandOpen ();
  G4String GetCurrentValue(G4UIcommand*);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandOpen (const G4VisCommandOpen&);
  G4VisCommandOpen& operator = (const G4VisCommandOpen&);
  G4UIcommand* fpCommand;
};

class G4VisCommandPlot: public G4VVisCommand {
public:
  G4VisCommandPlot ();
  virtual ~G4VisCommandPlot ();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandPlot (const G4VisCommandPlot&);
  G4VisCommandPlot& operator = (const G4VisCommandPlot&);
  G4UIcommand* fpCommand;
};

class G4VisCommandSpecify: public G4VVisCommand {
public:
  G4VisCommandSpecify ();
  virtual ~G4VisCommandSpecify ();
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandSpecify (const G4VisCommandSpecify&);
  G4VisCommandSpecify& operator = (const G4VisCommandSpecify&);
  G4UIcommand* fpCommand;
};

#endif
