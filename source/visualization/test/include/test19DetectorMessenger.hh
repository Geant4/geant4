// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test19DetectorMessenger.hh,v 1.2 1999-12-15 14:54:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Detector Construction Messenger for visualization testing.
// John Allison 25th April 1997

#ifndef test19DetectorMessenger_hh
#define test19DetectorMessenger_hh

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4VUserDetectorConstruction.hh"

class G4UIcommand;
class test19DetectorConstruction;

class test19DetectorMessenger: public G4UImessenger
{
public:
  test19DetectorMessenger (test19DetectorConstruction* test19Det);
  ~test19DetectorMessenger ();
  void SetNewValue (G4UIcommand* command, G4String newValues);
private:
  test19DetectorConstruction* test19Detector;
  G4UIcommand* fpTest19DetCommandDirectory;
  G4UIcommand* fpDetectorCommand;
};

#endif
