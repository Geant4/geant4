#ifndef TST50PROTONEEDL_HH
#define TST50PROTONEEDL_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class Tst50ProtonEEDL : public G4VPhysicsConstructor {

public: 

  Tst50ProtonEEDL(const G4String& name = "proton-eedl");
  
  virtual ~Tst50ProtonEEDL();
  
  // This method is dummy for physics
  virtual void ConstructParticle() {};
  
  virtual void ConstructProcess();
};

#endif

