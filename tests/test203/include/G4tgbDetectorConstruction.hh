// This class is the link between the HarpDD/Detrep packages and the G4RunManager
// It passes to the G4RunManager the top volume in the hierarchy constructed in G4tgbVolumeMgr:
#ifndef G4tgbDetectorConstruction_H
#define G4tgbDetectorConstruction_H 1
#include "globals.hh"

class G4VPhysicalVolume;

class G4tgbDetectorConstruction
{
 public:
  G4tgbDetectorConstruction();
  ~G4tgbDetectorConstruction();
  
 public:
  G4VPhysicalVolume* Construct();

};

#endif

