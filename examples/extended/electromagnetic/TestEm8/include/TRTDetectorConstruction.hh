
#ifndef TRTDetectorConstruction_h
#define TRTDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "Randomize.hh"

class TRTTrackerSD;
class SCTTrackerSD;
class HepJamesRandom;

class TRTDetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    TRTDetectorConstruction();
    ~TRTDetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

  private:
    TRTTrackerSD *TRTStrawSD;
    SCTTrackerSD *SCTStripSD;
    HepJamesRandom theJamesEngine;    
};

#endif
