#ifndef TST50PROTONSTANDARD_HH
#define TST50PROTONSTANDARD_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class Tst50ProtonStandard : public G4VPhysicsConstructor {

public: 

  Tst50ProtonStandard(const G4String& name = "proton-standard");
  
  virtual ~Tst50ProtonStandard();
  
  // This method is dummy for physics
  virtual void ConstructParticle() {};
  
  virtual void ConstructProcess();
};

#endif

