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
