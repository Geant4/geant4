///////////////////////////////////////////////////////////////////////////////
// File: HcalTestBeam96DetectorConstruction.hh
// Modifications: 23/08/00
// Constructs the 96 test beam configuration of HCAL
///////////////////////////////////////////////////////////////////////////////

#ifndef HcalTestBeam96DetectorConstruction_h
#define HcalTestBeam96DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class HcalTestBeam96DetectorConstruction : public G4VUserDetectorConstruction {
public:
  HcalTestBeam96DetectorConstruction();
  ~HcalTestBeam96DetectorConstruction();

public:
  G4VPhysicalVolume* Construct();
};

#endif
