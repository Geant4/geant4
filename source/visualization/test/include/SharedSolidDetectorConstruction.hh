#ifndef SharedSolidDetectorConstruction_h
#define SharedSolidDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"

class SharedSolidDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    SharedSolidDetectorConstruction();
    ~SharedSolidDetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();
};

#endif
