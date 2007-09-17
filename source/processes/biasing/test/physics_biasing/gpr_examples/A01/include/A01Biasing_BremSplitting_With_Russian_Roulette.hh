#ifndef A0BIASING_BREMSPLITTING_WITH_RUSSIAN_ROULETTE_HH
#define A0BIASING_BREMSPLITTING_WITH_RUSSIAN_ROULETTE_HH

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
#include "G4GammaConversion.hh"
#include "G4ComptonScattering.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4Gamma.hh"

using namespace G4GPRBuilder;

class A01Biasing_BremSplitting_With_Russian_Roulette : public G4VUserPhysicsBiasing {

public:
  
  void ConstructBiasing() 
  {
    // Simple case 
    G4GPRBiasingConfig electronCfg;
    electronCfg.SelectVProcess<G4eBremsstrahlung>();
    electronCfg.SelectParticle<G4Electron>();

    AddBiasing<G4GPRProcessLists::DiscreteDoIt>("Uniform Brem Splitting", 
						&A01BremSplittingFunctions::BremSplitting, electronCfg);

    G4GPRBiasingConfig gammaCfg;
    gammaCfg.SelectVProcess<G4GammaConversion>();
    gammaCfg.SelectVProcess<G4ComptonScattering>();
    gammaCfg.SelectVProcess<G4PhotoElectricEffect>();
    gammaCfg.SelectParticle<G4Gamma>();

    AddBiasing<G4GPRProcessLists::DiscreteDoIt>("Roulette", 
						&A01BremSplittingFunctions::Roulette, gammaCfg);

  }
};
#endif

