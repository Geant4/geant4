///////////////////////////////////////////////////////////////////////////////
// File: CCalDetectorConstruction.hh
// Description: Constructs the 96 test beam configuration of HCAL
///////////////////////////////////////////////////////////////////////////////

#ifndef CCalDetectorConstruction_h
#define CCalDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class CCalDetectorConstruction : public G4VUserDetectorConstruction {
public:
  CCalDetectorConstruction();
  ~CCalDetectorConstruction();

public:
  G4VPhysicalVolume* Construct();
};

#endif
