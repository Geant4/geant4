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
