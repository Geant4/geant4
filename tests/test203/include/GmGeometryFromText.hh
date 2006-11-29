#ifndef GmGeometryFromText_HH
#define GmGeometryFromText_HH

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

class GmGeometryFromText :public G4VUserDetectorConstruction{
 public:
  GmGeometryFromText();
  ~GmGeometryFromText(){};

  G4VPhysicalVolume* Construct();
};

#endif
