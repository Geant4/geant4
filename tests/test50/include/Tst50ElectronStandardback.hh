#ifndef TST50ELECTRONSTANDARDBACK_HH
#define TST50ELECTRONSTANDARDBACK_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class Tst50ElectronStandardback : public G4VPhysicsConstructor {

public: 

  Tst50ElectronStandardback(const G4String& name = "electron-standard-back");
  
  virtual ~Tst50ElectronStandardback();
  
  // This method is dummy for physics
  virtual void ConstructParticle() {};
  
  virtual void ConstructProcess();
};

#endif

