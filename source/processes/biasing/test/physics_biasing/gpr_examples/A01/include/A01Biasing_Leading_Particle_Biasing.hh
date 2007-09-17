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
#include "A01LeadingParticleBiasing_EM.hh"
#include "A01LeadingParticleBiasing_Hadronic.hh"
#include "G4GPRProcessWrappers.hh"
#include "G4GPRBiasingConfig.hh"
#include "G4HadronicProcess.hh"

using namespace G4GPRBuilder;
using namespace A01Triggers;

class A01Biasing_Leading_Particle_Biasing : public G4VUserPhysicsBiasing {

public:
  
  void ConstructBiasing() 
  {/*
    // Simple case 
    //AddBiasingWithTrigger<Particle, Process list, Trigger Type>(Name, Process Index, biasing function, trigger function)
    AddTriggeredBiasing<G4Gamma, G4GPRProcessLists::DiscreteDoIt, G4GPRTriggerTypes::Geometry::NewVolume>
      ("LeadingParticlebiasing", 1, &A01LeadingParticleBiasing_EM::SimpleEM_Conv, &CalorimeterTrigger);

    AddTriggeredBiasing<G4Electron, G4GPRProcessLists::DiscreteDoIt, G4GPRTriggerTypes::Geometry::NewVolume>
      ("LeadingParticlebiasing", 3, &A01LeadingParticleBiasing_EM::SimpleEM, &CalorimeterTrigger);

    AddTriggeredBiasing<G4Positron, G4GPRProcessLists::DiscreteDoIt, G4GPRTriggerTypes::Geometry::NewVolume>
      ("LeadingParticlebiasing", 3, &A01LeadingParticleBiasing_EM::SimpleEM, &CalorimeterTrigger);
   */
    // Slightly more complex configuration case - want to bias all G4HadronicProcess's
    G4GPRBiasingConfig hadronicConfig;
    hadronicConfig.SelectAllParticles();
    hadronicConfig.SelectVProcess<G4HadronicProcess>();
    
    AddTriggeredBiasing<G4GPRProcessLists::DiscreteDoIt, G4GPRTriggerTypes::Stepping::StartStep>
      ("LeadingParticlebiasing_Hadronic", 
       &A01LeadingParticleBiasing_Hadronic::Biasing, 
       &Hadronic_LeadingParticleBiasing_Trigger, hadronicConfig);
      
  }

};
#endif

