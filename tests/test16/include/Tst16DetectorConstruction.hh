// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst16DetectorConstruction.hh,v 1.1 1999-11-18 14:58:16 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst16DetectorConstruction_h
#define Tst16DetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class Tst16DetectorMessenger;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class Tst16DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    Tst16DetectorConstruction();
    ~Tst16DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();
     void SelectMaterial(G4String val);

  private:
     void SelectMaterialPointer();

     G4LogicalVolume*   simpleBoxLog;
     G4Material* Air;
     G4Material* Al;
     G4Material* Pb;
     G4Material* selectedMaterial;
     G4String materialChoice;
     Tst16DetectorMessenger * detectorMessenger;
};

#endif

