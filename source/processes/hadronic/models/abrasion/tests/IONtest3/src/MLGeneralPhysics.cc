////////////////////////////////////////////////////////////////////////////////
//
#include "MLGeneralPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>   

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
// Bosons
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
////////////////////////////////////////////////////////////////////////////////
//
MLGeneralPhysics::MLGeneralPhysics (const G4String& name)
  :G4VPhysicsConstructor(name)
{}
////////////////////////////////////////////////////////////////////////////////
//
MLGeneralPhysics::~MLGeneralPhysics ()
{}
////////////////////////////////////////////////////////////////////////////////
//
void MLGeneralPhysics::ConstructParticle ()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();  
}
////////////////////////////////////////////////////////////////////////////////
//

void MLGeneralPhysics::ConstructProcess()
{
  // Add Decay Process
  theParticleIterator->reset();
  while((*theParticleIterator)()){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (fDecayProcess.IsApplicable(*particle)) { 
      pmanager ->AddProcess(&fDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(&fDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(&fDecayProcess, idxAtRest);
    }
  }
}
////////////////////////////////////////////////////////////////////////////////

