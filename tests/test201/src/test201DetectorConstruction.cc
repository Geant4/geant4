// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test201DetectorConstruction.cc,v 1.2 1999-12-15 14:55:02 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Detector Construction for visualization testing.
// John Allison 24th April 1997

#include "test201DetectorConstruction.hh"

#include "test201DetectorMessenger.hh"
#include "G4UImanager.hh"

test201DetectorConstruction::test201DetectorConstruction ():
fpDetector (0) {
  new test201DetectorMessenger (this);
}

test201DetectorConstruction::~test201DetectorConstruction () {}

G4VPhysicalVolume* test201DetectorConstruction::Construct () {

  if (!fpDetector) {
    G4cout << "Detector not established - constructing default detector."
         << G4endl;
    G4UImanager* UI = G4UImanager::GetUIpointer ();
    UI -> ApplyCommand("/test201/detector 4");  // Sets fpDetector.
  }
  return fpDetector;
}
