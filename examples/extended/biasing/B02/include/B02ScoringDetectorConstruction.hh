#ifndef B02ScoringDetectorConstruction_hh 
#define B02ScoringDetectorConstruction_hh  B02ScoringDetectorConstruction_hh 
#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume; 
class B02ScoringDetectorConstruction : public G4VUserDetectorConstruction {
public:
  B02ScoringDetectorConstruction(){};
  ~B02ScoringDetectorConstruction(){};

  G4VPhysicalVolume* Construct();
  
};


#endif
