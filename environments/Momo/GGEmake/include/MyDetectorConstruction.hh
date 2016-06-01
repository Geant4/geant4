// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: MyDetectorConstruction.hh,v 1.2 1998/11/12 10:50:05 yhajime Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef MyDetectorConstruction_H
#define MyDetectorConstruction_H 1

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    MyDetectorConstruction();
    ~MyDetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

};

#endif

