#ifndef _DETECTORCONSTRUCTION_H_
#define _DETECTORCONSTRUCTION_H_

#include "G4VUserDetectorConstruction.hh"

class DetectorConstruction : public G4VUserDetectorConstruction {
private:
   G4VPhysicalVolume *World;
public:
 
   DetectorConstruction(G4VPhysicalVolume *setWorld = 0) {
   
      World = setWorld;
   }   

  G4VPhysicalVolume *Construct() {
  
     return World;
  }
};

#endif
