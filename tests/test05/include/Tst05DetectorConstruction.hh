
#ifndef Tst05DetectorConstruction_h
#define Tst05DetectorConstruction_h 1
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Timer.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Assembly.hh"
#include "G4AssemblyCreator.hh"
#include "globals.hh"

class Tst05DetectorMessenger;

class Tst05DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  Tst05DetectorConstruction();
  ~Tst05DetectorConstruction();
  
public:
  G4VPhysicalVolume* Construct();
  void SelectDetector(G4String val);
  
private:
  G4double expHall_x;
  G4double expHall_y;
  G4double expHall_z;
  
  G4double calBox_x;
  G4double calBox_y;
  G4double calBox_z;
  G4double rotAngle;
  G4double calPos;

  G4double trackerRadius;
  G4double trackerHight;
  G4double trackerPos;

  G4int detectorChoice;
  Tst05DetectorMessenger* detectorMessenger;
};

#endif

