#ifndef A0BIASING_LEADING_PARTICLE_BIASING_HH
#define A0BIASING_LEADING_PARTICLE_BIASING_HH

#include "A01Triggers.hh"

#include "G4VUserPhysicsBiasing.hh"
#include "G4GPRBuilder.hh"
#include "G4GPRTriggerTypes.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Transportation.hh"
#include "A01LeadingParticleBiasing.hh"
#include "G4GPRProcessWrappers.hh"

using namespace G4GPRBuilder;
using namespace A01Triggers;

class A01Biasing_Leading_Particle_Biasing : public G4VUserPhysicsBiasing {

public:
  
  void ConstructBiasing() 
  {
    // Simple case 
    //AddRelayWithTrigger<Particle, Process list, Trigger Type>(Name, Process Index, biasing function, trigger function)
    AddRelayWithTrigger<G4Gamma, G4GPRProcessLists::DiscreteDoIt, G4GPRTriggerTypes::Geometry::NewVolume>
      ("LeadingParticlebiasing", 1, &A01LeadingParticleBiasing::SimpleEM_Conv, &CalorimeterTrigger);

    AddRelayWithTrigger<G4Electron, G4GPRProcessLists::DiscreteDoIt, G4GPRTriggerTypes::Geometry::NewVolume>
      ("LeadingParticlebiasing", 3, &A01LeadingParticleBiasing::SimpleEM, &CalorimeterTrigger);

    AddRelayWithTrigger<G4Positron, G4GPRProcessLists::DiscreteDoIt, G4GPRTriggerTypes::Geometry::NewVolume>
      ("LeadingParticlebiasing", 3, &A01LeadingParticleBiasing::SimpleEM, &CalorimeterTrigger);
  }
};
#endif

