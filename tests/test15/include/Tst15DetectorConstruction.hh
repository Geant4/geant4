// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst15DetectorConstruction.hh,v 1.3 2000-02-23 10:50:17 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst15DetectorConstruction_h
#define Tst15DetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class Tst15DetectorMessenger;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class Tst15DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    Tst15DetectorConstruction();
    ~Tst15DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();
     void SelectMaterial(G4String val);

  private:
     void SelectMaterialPointer();

     G4LogicalVolume*   simpleBoxLog;
     G4Material* Air;
     G4Material* Al;
     G4Material* Pb;
     G4Material* U;
     G4Material* selectedMaterial;
     G4String materialChoice;
     Tst15DetectorMessenger * detectorMessenger;
};

#endif

