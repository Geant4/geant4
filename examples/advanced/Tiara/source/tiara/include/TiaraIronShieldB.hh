// $Id: TiaraIronShieldB.hh,v 1.2 2003-06-16 17:06:46 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
//
// Class TiaraIronShieldB
//

#ifndef TiaraIronShieldB_hh
#define TiaraIronShieldB_hh TiaraIronShieldB_hh

#include "TiaraVComponent.hh"


class TiaraDimensions;
class TiaraMaterials;

class TiaraIronShieldB : public TiaraVComponent {
public:
  TiaraIronShieldB(TiaraMaterials &mfac, 
		   const TiaraDimensions &tiaraDimensions);
  ~TiaraIronShieldB();

  virtual TiaraParts GetParts();

private:
  TiaraParts fTiaraParts;
};

#endif
