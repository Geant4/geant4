// $Id: TiaraVComponent.hh,v 1.3 2003-06-18 16:40:25 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
//
// Class TiaraVComponent
//

#ifndef TiaraVComponent_hh
#define TiaraVComponent_hh TiaraVComponent_hh

#include <vector>
#include "G4ThreeVector.hh"
#include "TiaraPart.hh"

class G4LogicalVolume;

typedef std::vector<TiaraPart> TiaraParts;

class TiaraVComponent {
public:
  TiaraVComponent();
  virtual ~TiaraVComponent();

  virtual TiaraParts GetParts() = 0;
};

#endif
