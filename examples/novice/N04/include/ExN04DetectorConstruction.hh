
#ifndef ExN04DetectorConstruction_h
#define ExN04DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4VPhysicalVolume;

class ExN04DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ExN04DetectorConstruction();
    ~ExN04DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();

  private:

#include "ExN04DetectorParameterDef.hh"

};

#endif

