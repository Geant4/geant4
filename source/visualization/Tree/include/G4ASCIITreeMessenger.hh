// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ASCIITreeMessenger.hh,v 1.3 2001-06-15 07:22:55 stanaka Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  5th April 2001
// A scene handler to dump geometry hierarchy to standard output as
//   ASCII stream.
// Based on a provisional G4ASCIITreeGraphicsScene (was in modeling).

#ifndef G4ASCIITREEMESSENGER_HH
#define G4ASCIITREEMESSENGER_HH

#include "G4UImessenger.hh"

class G4UIcommand;
class G4UIcmdWithAnInteger;
class G4ASCIITree;

class G4ASCIITreeMessenger: public G4UImessenger {
public:
  G4ASCIITreeMessenger(G4ASCIITree*);
  virtual ~G4ASCIITreeMessenger();
  G4String GetCurrentValue (G4UIcommand* command);
  void SetNewValue (G4UIcommand* command, G4String newValue);
private:
  G4ASCIITree* fpASCIITree;
  G4UIcommand* fpDirectory;
  G4UIcmdWithAnInteger* fpCommandVerbose;
};

#endif
