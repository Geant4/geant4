// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4VisToOldVisCommands.hh,v 1.5 2001-02-23 15:43:20 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Implements some /vis/ commands as /vis~/ temporarily.
// John Allison  9th December 1998.

#ifndef G4VISTOOLDVISCOMMANDS_HH
#define G4VISTOOLDVISCOMMANDS_HH

#include "G4UImessenger.hh"
#include "globals.hh"

class G4VisToOldVisCommands: public G4UImessenger {
public:
  G4VisToOldVisCommands ();
  ~G4VisToOldVisCommands ();
  void SetNewValue (G4UIcommand* command, G4String newValues);
};

#endif
