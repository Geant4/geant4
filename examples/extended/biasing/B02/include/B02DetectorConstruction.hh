
#ifndef B02DetectorConstruction_hh
#define B02DetectorConstruction_hh B02DetectorConstruction_hh

class G4VPhysicalVolume;

#include "G4VUserDetectorConstruction.hh"

class B02DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  B02DetectorConstruction();
  ~B02DetectorConstruction();
  
  G4VPhysicalVolume* Construct();
  G4VPhysicalVolume* GetWorldVolume(){return fWorldVolume;}
private:
  G4VPhysicalVolume* fWorldVolume;
};

#endif

