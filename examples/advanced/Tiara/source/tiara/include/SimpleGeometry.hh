#ifndef SimpleGeometry_hh 
#define SimpleGeometry_hh SimpleGeometry_hh

#include "G4VUserDetectorConstruction.hh"

class SimpleGeometry : public G4VUserDetectorConstruction {
public:
  SimpleGeometry(G4VPhysicalVolume *); 
  ~SimpleGeometry();
  G4VPhysicalVolume* Construct();
private:
  G4VPhysicalVolume *fVolume;
};


#endif
