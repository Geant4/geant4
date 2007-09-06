#ifndef A01BIASING_HH
#define A01BIASING_HH

#include "G4VUserPhysicsBiasing.hh"
//#include "G4GPRBuilder.hh"

//using namespace G4GPRBuilder;
/*
G4bool primaryElectronTrigger(const G4Track& track)
{
  return track->GetDef() == G4Electron::Definition();
}

G4bool calorimeterTrigger(const G4Track& track, const G4Step& step)
{
  return (track.GetVolume() == "cellPhysical");
}

*/
class A01Biasing : public G4VUserPhysicsBiasing {

public:

  void ConstructBiasing();
    /*
  {

    // Override default A01 physics list for positrons. Just for fun, make the positrons use
    // the low energy physics list by default
    G4String newPhysicsList("New Default");
    CreateDefaultPhysicsList<G4Positron>(newPhysicsList);

    AddProcess<G4Positron>(new G4Transportation, x, newPhysicsList);
    AddProcess<G4Positron>(new G4MultipleScattering, -1, 1, 1, newPhysicsList);
    AddProcess<G4Positron>(new G4eIonisation,        -1, 2, 2, newPhysicsList);
    AddProcess<G4Positron>(new G4eBremsstrahlung,    -1, 3, 3, newPhysicsList);
    AddProcess<G4Positron>(new G4eplusAnnihilation,   0,-1, 4, newPhysicsList);

    // Create new physics list for gammas which is only active when in calorimeter volume
    CreatePhysicsListWithTrigger<G4Gamma, G4GPRTriggering::Geometry::StartBoundary>("Calorimeter Physics List", &calorimeterTrigger);
    G4GPRProcessConfig cfg1;
    cfg.SelectParticle<G4Gamma>();
    cfg.SelectIdx(G4GPRProcessPlacement::First);

    AddProcess<G4GPRProcessLists::Discrete>(new G4Transportation, cfg);

    G4GPRProcessConfig cfg2;
    cfg2.SelectParticle<G4Gamma>();
    cfg2.SelectIdx(G4GPRProcessPlacement::Last);
    AddProcess<G4GPRProcessLists::Discrete>(myGammaProcess, cfg2);

    // Uniform bremsstrahlung splitting for primary electrons
    G4GPRSingleRelayConfig<G4GPRProcessTypes::DiscreteDoIt> relayConfig;
    relayConfig.AddParticle<G4Electron>();
    relayConfig.AddtoProcess(3);

    // Further development - could select processes to attach to based on template parameter
    // relayConfig.AddToProcess<G4eBremsstrahlung>();

    AddSingleRelayWithTrigger<G4GPRTriggering::Geometry::StartTracking>(new MyBremSplitting, relayConfig, &primaryTrackTrigger);
    
    // Hadronic leading particle biasing for hadrons with energy less than 5GeV
    G4GPRSingleRelayConfig leadParticleBiasingConfig;
    leadParticleBiasingConfig.AddParticles(&myLeadParticleBiasingSelection);
    leadParticleBiasing.AddToProcess<G4Lala>();
    
    G4GPRMinEnergyTrigger* trigger = new G4GPRMinEnergyTrigger;
    trigger->SetMinimumEnergy(5*GeV);
    AddSingleRelayWithTrigger(&NewHadLeadBiasing, leadParticleBiasingConfig, trigger);
    
    //future development - combo trigger, less than 5gev and lala

    // Attempt at implicit capture
    G4GPRMultiRelayConfig implicitCaptureConfig;
    implicitCaptureConfig.AddParticle<G4Neutron>();
    implicitCaptureConfig.AddProcess<>;
    implicitCaptureConfig.AddProcess<>;
    implicitCaptureConfig.AddProcess<>;
    
    G4GPRImplicitCapture* implicitCapture = new G4GPRImplicitCapture;
    AddMultiRelay(implicitCapture, implicitCaptureConfig);
  }
    */

};
#endif

