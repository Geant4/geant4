#ifndef PCTBinaryCascadeWrapper_hh
#define PCTBinaryCascadeWrapper_hh

#include "G4BinaryCascade.hh"
#include "G4VPreCompoundModel.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "PCTProjectile.hh"

class PCTBinaryCascadeWrapper
{
public:
  PCTBinaryCascadeWrapper() {}
  ~PCTBinaryCascadeWrapper() {}

  void SetDeExcitation(G4VPreCompoundModel * model)
  {
    theCascade.SetDeExcitation(model);
    return;
  }

  G4ReactionProductVector * DeExcite(const PCTProjectile * theProjectile,
				     const G4int A, const G4int Z);

  

private:

  G4BinaryCascade theCascade;
};
#endif
