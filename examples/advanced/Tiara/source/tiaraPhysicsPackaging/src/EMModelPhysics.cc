#include "EMModelPhysics.hh"
#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"


EMModelPhysics::
EMModelPhysics(const G4String& name) : G4VPhysicsConstructor(name) {}

EMModelPhysics::
~EMModelPhysics() {}

void EMModelPhysics::
ConstructParticle()
{
  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
}

void EMModelPhysics::
ConstructProcess()
{
  theEMModelPhysics.Build();
}



// 2002 by J.P. Wellisch
