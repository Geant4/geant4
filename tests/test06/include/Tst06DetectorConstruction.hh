// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst06DetectorConstruction.hh,v 1.2 1999-12-15 14:54:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef Tst06DetectorConstruction_h
#define Tst06DetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class Tst06DetectorMessenger;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class Tst06DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    Tst06DetectorConstruction();
    ~Tst06DetectorConstruction();

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
     Tst06DetectorMessenger * detectorMessenger;
};

#endif

