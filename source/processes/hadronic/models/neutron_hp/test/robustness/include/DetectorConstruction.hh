// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: DetectorConstruction.hh,v 1.1 2001-05-29 19:17:32 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction();

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
     DetectorMessenger * detectorMessenger;
};

#endif

