// $Id: TiaraVComponent.hh,v 1.2 2003-06-16 17:06:46 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
//
// Class TiaraVComponent
//

#ifndef TiaraVComponent_hh
#define TiaraVComponent_hh TiaraVComponent_hh

#include "g4std/vector"
#include "G4ThreeVector.hh"
#include "TiaraPart.hh"

class G4LogicalVolume;

typedef G4std::vector<TiaraPart> TiaraParts;

class TiaraVComponent {
public:
  TiaraVComponent();
  virtual ~TiaraVComponent();

  virtual TiaraParts GetParts() = 0;
};

#endif
