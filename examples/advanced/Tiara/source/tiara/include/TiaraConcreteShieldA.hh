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
