#include "EMPhysics.hh"
#include "globals.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"


EMPhysics::
EMPhysics(const G4String& name) : G4VPhysicsConstructor(name) {}

EMPhysics::
~EMPhysics() {}

void EMPhysics::
ConstructParticle()
{
  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
}

void EMPhysics::
ConstructProcess()
{
  theEMPhysics.Build();
}



// 2002 by J.P. Wellisch
