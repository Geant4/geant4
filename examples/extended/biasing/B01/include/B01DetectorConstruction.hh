
#ifndef B01DetectorConstruction_hh
#define B01DetectorConstruction_hh B01DetectorConstruction_hh

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class B01DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  B01DetectorConstruction();
  ~B01DetectorConstruction();
  
  G4VPhysicalVolume* Construct();
  G4VPhysicalVolume* GetWorldVolume(){return fWorldVolume;}
private:
  G4VPhysicalVolume* fWorldVolume;
};

#endif

