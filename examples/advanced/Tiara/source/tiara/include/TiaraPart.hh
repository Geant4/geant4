// $Id: TiaraPart.hh,v 1.2 2003-06-16 17:06:46 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
//
// Class TiaraPart
//

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
