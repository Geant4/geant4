// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyDetectorMessenger.hh,v 1.1 1999-04-16 10:32:29 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyDetectorMessenger_h
#define MyDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4UIcommand.hh"

class MyDetectorConstruction;

class MyDetectorMessenger: public G4UImessenger
{
public:
  MyDetectorMessenger(MyDetectorConstruction * myDet);
  ~MyDetectorMessenger ();
  void SetNewValue(G4UIcommand * command,G4String newValues);
private:
  MyDetectorConstruction * myDetector;
  G4UIcommand* fpMyDetCommandDirectory;
  G4UIcommand* fpCalMaterialCommand;
};

#endif

