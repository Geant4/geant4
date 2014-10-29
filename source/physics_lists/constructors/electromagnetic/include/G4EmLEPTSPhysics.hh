#ifndef __G4EmLEPTSPhysics__
#define __G4EmLEPTSPhysics__

#include "G4VPhysicsConstructor.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
// Geant461: #include "g4std/iomanip"
#include "G4Decay.hh"
#include "G4UserSpecialCuts.hh"
#include "G4StepLimiter.hh"
#include "globals.hh"


class G4EmLEPTSPhysics : public G4VPhysicsConstructor
{
 public:
  G4EmLEPTSPhysics(const G4String& name="G4EmLEPTSPhysics");
  virtual ~G4EmLEPTSPhysics(){};

  virtual void ConstructParticle();
  virtual void ConstructProcess();

 protected:
  // Physics

};


#endif
