// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyDetectorConstruction.hh,v 1.1 1999-01-07 16:04:59 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef MyDetectorConstruction_h
#define MyDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    MyDetectorConstruction();
    ~MyDetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

};

#endif

