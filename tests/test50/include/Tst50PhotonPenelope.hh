
#ifndef TST50PHOTONPENELOPE_HH
#define TST50PHOTONPENELOPE_HH 1

#include "G4VPhysicsConstructor.hh"
#include "globals.hh"

class Tst50PhotonPenelope : public G4VPhysicsConstructor {

public: 

  Tst50PhotonPenelope(const G4String& name = "photon-penelope");
  
  virtual ~Tst50PhotonPenelope();
  
  // This method is dummy for physics
  virtual void ConstructParticle() {};
  
  virtual void ConstructProcess();
};

#endif

