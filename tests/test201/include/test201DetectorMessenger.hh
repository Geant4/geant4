// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test201DetectorMessenger.hh,v 1.1 1999-01-08 16:35:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Detector Construction Messenger for visualization testing.
// John Allison 25th April 1997

#ifndef test201DetectorMessenger_hh
#define test201DetectorMessenger_hh

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4VUserDetectorConstruction.hh"

class G4UIcommand;
class test201DetectorConstruction;

class test201DetectorMessenger: public G4UImessenger
{
public:
  test201DetectorMessenger (test201DetectorConstruction* test201Det);
  ~test201DetectorMessenger ();
  void SetNewValue (G4UIcommand* command, G4String newValues);
private:
  test201DetectorConstruction* test201Detector;
  G4UIcommand* fpTest201DetCommandDirectory;
  G4UIcommand* fpDetectorCommand;
};

#endif
