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
