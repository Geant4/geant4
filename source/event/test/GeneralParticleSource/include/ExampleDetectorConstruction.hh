
#ifndef ExampleDetectorConstruction_h
#define ExampleDetectorConstruction_h 1
#include "G4PVPlacement.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Assembly.hh"
#include "G4AssemblyCreator.hh"

class ExampleDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ExampleDetectorConstruction();
    ~ExampleDetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

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
};

#endif

