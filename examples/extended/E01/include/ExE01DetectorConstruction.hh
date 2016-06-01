
#ifndef ExE01DetectorConstruction_h
#define ExE01DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class ExE01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ExE01DetectorConstruction();
    ~ExE01DetectorConstruction();

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

