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
// $Id: G4VisCommandsCompound.hh,v 1.7 2001-07-11 10:09:15 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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

class G4VisCommandDrawVolume: public G4VVisCommand {
public:
  G4VisCommandDrawVolume ();
  virtual ~G4VisCommandDrawVolume ();
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandDrawVolume (const G4VisCommandDrawVolume&);
  G4VisCommandDrawVolume& operator = (const G4VisCommandDrawVolume&);
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandOpen: public G4VVisCommand {
public:
  G4VisCommandOpen ();
  virtual ~G4VisCommandOpen ();
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4VisCommandOpen (const G4VisCommandOpen&);
  G4VisCommandOpen& operator = (const G4VisCommandOpen&);
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
  G4UIcmdWithAString* fpCommand;
};

#endif
