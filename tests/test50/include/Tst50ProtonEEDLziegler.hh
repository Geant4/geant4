#ifndef TST50PROTONEEDLziegler_HH
#define TST50PROTONEEDLziegler_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class Tst50ProtonEEDLziegler : public G4VPhysicsConstructor {

public: 

  Tst50ProtonEEDLziegler(const G4String& name = "proton-eedl-ziegler");
  
  virtual ~Tst50ProtonEEDLziegler();
  
  // This method is dummy for physics
  virtual void ConstructParticle() {};
  
  virtual void ConstructProcess();
};

#endif

