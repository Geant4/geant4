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
