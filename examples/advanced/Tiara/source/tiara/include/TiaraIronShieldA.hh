// $Id: TiaraIronShieldA.hh,v 1.2 2003-06-16 17:06:46 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
//
// Class TiaraIronShieldA
//

#ifndef TiaraIronShieldA_hh
#define TiaraIronShieldA_hh TiaraIronShieldA_hh

#include "TiaraVComponent.hh"

class TiaraDimensions;
class TiaraMaterials;

class TiaraIronShieldA : public TiaraVComponent {
public:
  TiaraIronShieldA(TiaraMaterials &mfac,
		   const TiaraDimensions &tiaraDimensions);
  ~TiaraIronShieldA();

  virtual TiaraParts GetParts();

private:
  TiaraParts fTiaraParts;
};

#endif
