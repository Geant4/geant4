#ifndef Tst28DetectorConstruction_h
#define Tst28DetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class Tst28DetectorMessenger;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class Tst28DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    Tst28DetectorConstruction();
    ~Tst28DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();
     void SelectMaterial(G4String val);

  private:
     void SelectMaterialPointer();

     G4LogicalVolume*   simpleBoxLog;
     G4Material* theH;
     G4Material* theSi;
     G4Material* theCu;
     G4Material* thePb;
     G4Material* theU;
     G4Material* theH2O;
     
     G4Material* selectedMaterial;
     G4String materialChoice;
     Tst28DetectorMessenger * detectorMessenger;
};

#endif

