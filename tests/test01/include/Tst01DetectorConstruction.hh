
#ifndef Tst01DetectorConstruction_h
#define Tst01DetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class Tst01DetectorMessenger;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class Tst01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    Tst01DetectorConstruction();
    ~Tst01DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();
     void SwitchDetector();
     void SelectDetector(G4String val);
     void SelectMaterial(G4String val);

  private:
     G4VPhysicalVolume* SelectDetector();
     void ConstructDetectors();
     void SelectMaterialPointer();

     G4LogicalVolume*   simpleBoxLog;
     G4VPhysicalVolume* simpleBoxDetector;
     G4VPhysicalVolume* honeycombDetector;
     G4Material* Air;
     G4Material* Al;
     G4Material* Pb;
     G4Material* selectedMaterial;
     G4int detectorChoice;
     G4String materialChoice;
     Tst01DetectorMessenger * detectorMessenger;
};

#endif

