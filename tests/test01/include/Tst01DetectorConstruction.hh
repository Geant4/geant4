
#ifndef Tst01DetectorConstruction_h
#define Tst01DetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4VSolid ;
class G4Material;
class Tst01DetectorMessenger;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class Tst01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

     Tst01DetectorConstruction();
    ~Tst01DetectorConstruction();

     G4VPhysicalVolume* Construct();

     void SelectDetector(G4String val) ;
     void SwitchDetector() ;

     void SelectMaterial(G4String val) ;

  // Select/Switch CSG/Boolean

     void SelectCSG(G4String name) ;
     void SwitchCSG() ;

     void SelectBoolean(G4String name) ;
     void SwitchBoolean() ;

  private:

  // Private functions

     G4VPhysicalVolume* SelectDetector();
     void ConstructDetectors();
     void SelectMaterialPointer();

  // Class members

     Tst01DetectorMessenger* detectorMessenger ;

     G4LogicalVolume*   simpleBoxLog ;
     G4VPhysicalVolume* simpleBoxDetector ;
     G4VPhysicalVolume* honeycombDetector ;
     G4VPhysicalVolume* fWorldPhysVol ;

     G4VSolid*          fTestCSG ;
     G4LogicalVolume*   fTestLog ;
     G4VPhysicalVolume* fTestVol ;

  // Daughter CSGs /logic/physics volumes

     G4VSolid*          fTestD1CSG ;
     G4LogicalVolume*   fTestD1Log ;
     G4VPhysicalVolume* fTestD1Vol ;

     G4VSolid*          fTestD2CSG ;
     G4LogicalVolume*   fTestD2Log ;
     G4VPhysicalVolume* fTestD2Vol ;

     G4Material* Air ;
     G4Material* Al ;
     G4Material* Pb ;
     G4Material* selectedMaterial ;

     G4int       detectorChoice ;
     G4String    materialChoice ;

     G4int       fChoiceCSG ;
     G4int       fChoiceBool ;
};

#endif

