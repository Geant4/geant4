#ifndef A0BIASING_BREMSPLITTING_DAUGHTEROFPRIMARY_ONLY_HH
#define A0BIASING_BREMSPLITTING_DAUGHTEROFPRIMARY_ONLY_HH

#include "A01Triggers.hh"

#include "G4VUserPhysicsBiasing.hh"
#include "G4GPRBuilder.hh"
#include "G4GPRTriggerTypes.hh"
#include "G4Positron.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4GPRBiasingConfig.hh"
#include "G4eBremsstrahlung.hh"
#include "A01BremSplittingFunctions.hh"

using namespace G4GPRBuilder;

class A01Biasing_BremSplitting_DaugherOfPrimary_Only : public G4VUserPhysicsBiasing {

public:
  
  void ConstructBiasing() 
  {
    // Simple case 
    G4GPRBiasingConfig electronCfg;
    electronCfg.SelectVProcess<G4eBremsstrahlung>();
    electronCfg.SelectParticle<G4Electron>();

    AddTriggeredBiasing<G4GPRProcessLists::DiscreteDoIt, G4GPRTriggerTypes::Tracking::StartTracking>
      ("Splitting For Primary Electron", &A01BremSplittingFunctions::BremSplitting, &A01Triggers::DaughterOfPrimaryTrigger, electronCfg);

  }
};
#endif

