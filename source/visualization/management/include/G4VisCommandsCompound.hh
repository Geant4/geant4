// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisCommandsCompound.hh,v 1.2 2000-05-15 11:49:02 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

// Compound /vis/ commands - John Allison  15th May 2000

#ifndef G4VISCOMMANDSCOMPOUND_HH
#define G4VISCOMMANDSCOMPOUND_HH

#include "G4VVisCommand.hh"

class G4VisCommandDrawVolume: public G4VVisCommand {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandDrawVolume ();
  ~G4VisCommandDrawVolume ();
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandOpen: public G4VVisCommand {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandOpen ();
  ~G4VisCommandOpen ();
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithAString* fpCommand;
};

class G4VisCommandSpecify: public G4VVisCommand {
public:
  // Uses compiler defaults for copy constructor and assignment.
  G4VisCommandSpecify ();
  ~G4VisCommandSpecify ();
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4UIcmdWithAString* fpCommand;
};

#endif
