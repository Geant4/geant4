// $Id: TiaraConcreteShieldB.hh,v 1.2 2003-06-16 17:06:45 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
//
// Class TiaraConcreteShieldB
//

#ifndef TiaraConcreteShieldB_hh
#define TiaraConcreteShieldB_hh TiaraConcreteShieldB_hh

#include "TiaraVComponent.hh"

class TiaraDimensions;
class TiaraMaterials;

class TiaraConcreteShieldB : public TiaraVComponent {
public:
  TiaraConcreteShieldB(TiaraMaterials &mfac,
		       const TiaraDimensions &tiaraDimensions);
  ~TiaraConcreteShieldB();

  virtual TiaraParts GetParts();

private:
  TiaraParts fTiaraParts;
};

#endif
