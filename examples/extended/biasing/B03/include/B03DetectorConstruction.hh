#ifndef B03DetectorConstruction_hh
#define B03DetectorConstruction_hh B03DetectorConstruction_hh
#include "globals.hh"

class G4VPhysicalVolume;
class G4VIStore;

#include "G4VUserDetectorConstruction.hh"

class B03DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  B03DetectorConstruction();
  ~B03DetectorConstruction();
  
  G4VPhysicalVolume* Construct();
  G4VIStore *GetIStore(){
    if (!fIStore) G4Exception("B03DetectorConstruction::fIStore empty!");
    return fIStore;
  }
private:
  G4VIStore *fIStore;
};

#endif

