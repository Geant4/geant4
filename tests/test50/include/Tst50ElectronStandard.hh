#ifndef TST50ELECTRONSTANDARD_HH
#define TST50ELECTRONSTANDARD_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class Tst50ElectronStandard : public G4VPhysicsConstructor {

public: 

  Tst50ElectronStandard(const G4String& name = "electron-standard");
  
  virtual ~Tst50ElectronStandard();
  
  // This method is dummy for physics
  virtual void ConstructParticle() {};
  
  virtual void ConstructProcess();
};

#endif

