// $Id: TiaraConcreteShieldA.hh,v 1.2 2003-06-16 17:06:45 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
//
// Class TiaraConcreteShieldA
//

#ifndef TiaraConcreteShieldA_hh
#define TiaraConcreteShieldA_hh TiaraConcreteShieldA_hh

#include "TiaraVComponent.hh"


class TiaraDimensions;
class TiaraMaterials;


class TiaraConcreteShieldA : public TiaraVComponent {
public:
  TiaraConcreteShieldA(TiaraMaterials &mfac,
		       const TiaraDimensions &tiaraDimensions);
  ~TiaraConcreteShieldA();

  virtual TiaraParts GetParts();

private:
  TiaraParts fTiaraParts;
};

#endif
