#ifndef TiaraPart_hh
#define TiaraPart_hh TiaraPart_hh

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

class G4LogicalVolume;

struct TiaraPart{
  TiaraPart() :
    logVol(0),
    pos(0,0,0),
    rot(0)
  {}
  G4LogicalVolume *logVol;
  G4ThreeVector pos;
  G4RotationMatrix *rot;
};

#endif 
