// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test19DetectorConstruction.cc,v 1.1 1999-04-16 10:32:39 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Detector Construction for visualization testing.
// John Allison 24th April 1997

#include "test19DetectorConstruction.hh"

#include "test19DetectorMessenger.hh"
#include "G4UImanager.hh"

test19DetectorConstruction::test19DetectorConstruction ():
fpDetector (0) {
  new test19DetectorMessenger (this);
}

test19DetectorConstruction::~test19DetectorConstruction () {}

G4VPhysicalVolume* test19DetectorConstruction::Construct () {

  if (!fpDetector) {
    G4cout << "Detector not established - constructing default detector."
         << endl;
    G4UImanager* UI = G4UImanager::GetUIpointer ();
    UI -> ApplyCommand("/test19det/detector 4");  // Sets fpDetector.
  }
  return fpDetector;
}
