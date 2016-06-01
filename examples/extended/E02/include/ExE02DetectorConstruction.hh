
#ifndef ExE02DetectorConstruction_h
#define ExE02DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class ExE02DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ExE02DetectorConstruction();
    ~ExE02DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

};

#endif

