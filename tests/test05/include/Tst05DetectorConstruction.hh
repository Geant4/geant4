// $Id: Tst05DetectorConstruction.hh,v 1.3 2000-02-25 16:56:41 gcosmo Exp $
// ------------------------------------------------------------

#ifndef Tst05DetectorConstruction_h
#define Tst05DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

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

