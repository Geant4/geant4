#include "G4HadronQEDBuilder.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "g4std/iomanip"   
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ProcessManager.hh"

#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4SigmaMinus.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4SigmaPlus.hh"
#include "G4AntiSigmaPlus.hh"
#include "G4XiMinus.hh"
#include "G4AntiXiMinus.hh"
#include "G4OmegaMinus.hh"
#include "G4AntiOmegaMinus.hh"

G4HadronQEDBuilder::G4HadronQEDBuilder() {}
G4HadronQEDBuilder::~G4HadronQEDBuilder(){}

void G4HadronQEDBuilder::
RegisterOne(G4ProcessManager* aP, G4MultipleScattering * aM, G4hIonisation* aI)
{
  aP->AddProcess(aI, ordInActive,2, 2);
  aP->AddProcess(aM);
  aP->SetProcessOrdering(aM, idxAlongStep, 1);
  aP->SetProcessOrdering(aM, idxPostStep, 1);
}

void G4HadronQEDBuilder::Build()
{
  // PionPlus
  RegisterOne(G4PionPlus::PionPlus()->GetProcessManager(), &thePionPlusMult, &thePionPlusIonisation);

  // PionMinus
  RegisterOne(G4PionMinus::PionMinus()->GetProcessManager(), &thePionMinusMult, &thePionMinusIonisation);

  // KaonPlus
  RegisterOne(G4KaonPlus::KaonPlus()->GetProcessManager(), &theKaonPlusMult, &theKaonPlusIonisation);

  // KaonMinus
  RegisterOne(G4KaonMinus::KaonMinus()->GetProcessManager(), &theKaonMinusMult, &theKaonMinusIonisation);

  // Proton
  RegisterOne(G4Proton::Proton()->GetProcessManager(), &theProtonMult, &theProtonIonisation);

  // anti-Proton
  RegisterOne(G4AntiProton::AntiProton()->GetProcessManager(), &theAntiProtonMult, &theAntiProtonIonisation);
    
  // SigmaMinus
  RegisterOne(G4SigmaMinus::SigmaMinus()->GetProcessManager(), &theSigmaMinusMult, &theSigmaMinusIonisation);

  // anti-SigmaMinus
  RegisterOne(G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager(), &theAntiSigmaMinusMult, &theAntiSigmaMinusIonisation);

  // SigmaPlus
  RegisterOne(G4SigmaPlus::SigmaPlus()->GetProcessManager(), &theSigmaPlusMult, &theSigmaPlusIonisation);

  // anti-SigmaPlus
  RegisterOne(G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager(), &theAntiSigmaPlusMult, &theAntiSigmaPlusIonisation);

  // XiMinus
  RegisterOne(G4XiMinus::XiMinus()->GetProcessManager(), &theXiMinusMult, &theXiMinusIonisation);

  // anti-XiMinus
  RegisterOne(G4AntiXiMinus::AntiXiMinus()->GetProcessManager(), &theAntiXiMinusMult, &theAntiXiMinusIonisation);

  // OmegaMinus
  RegisterOne(G4OmegaMinus::OmegaMinus()->GetProcessManager(), &theOmegaMinusMult, &theOmegaMinusIonisation);

  // anti-OmegaMinus
  RegisterOne(G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager(), &theAntiOmegaMinusMult, &theAntiOmegaMinusIonisation);
}





// 2002 by J.P. Wellisch
