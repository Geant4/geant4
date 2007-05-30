//******************************************************************************
// DetectorConstruction.hh
//
// 1.00 RMK, LANL, MAR-2002:  First version.
//******************************************************************************
//
#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction();
    ~DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

};

#endif
