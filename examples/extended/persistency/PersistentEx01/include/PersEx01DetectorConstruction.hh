
#ifndef PersEx01DetectorConstruction_h
#define PersEx01DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class PersEx01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    PersEx01DetectorConstruction();
    ~PersEx01DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

};

#endif

