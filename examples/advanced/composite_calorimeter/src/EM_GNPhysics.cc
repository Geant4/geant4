#include "EM_GNPhysics.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"


EM_GNPhysics::
EM_GNPhysics(const G4String& name) : G4VPhysicsConstructor(name) {}

EM_GNPhysics::
~EM_GNPhysics() {}

void EM_GNPhysics::
ConstructParticle()
{
  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
}

void EM_GNPhysics::
ConstructProcess()
{
  theEMPhysics.Build();
  #ifndef NO_ELECTRO_AND_GAMMA_NUCLEAR
  theGNPhysics.Build();
  #endif
}



// 2002 by J.P. Wellisch
