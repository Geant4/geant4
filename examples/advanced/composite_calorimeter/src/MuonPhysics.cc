#include "MuonPhysics.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4TauMinus.hh"
#include "G4TauPlus.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoTau.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"

#include "G4ProcessManager.hh"

#include "globals.hh"
#include "G4ios.hh"

MuonPhysics::MuonPhysics(const G4String& name)
                   :  G4VPhysicsConstructor(name)
{
}

MuonPhysics::~MuonPhysics()
{
  if(wasActivated)
  {
  G4ProcessManager * pManager = 0;
  pManager = G4MuonPlus::MuonPlus()->GetProcessManager();
  if(pManager) pManager->RemoveProcess(&fMuPlusIonisation);
  if(pManager) pManager->RemoveProcess(&fMuPlusBremsstrahlung);
  if(pManager) pManager->RemoveProcess(&fMuPlusPairProduction);
  if(pManager) pManager->RemoveProcess(&fMuPlusMultipleScattering);
  pManager = G4MuonMinus::MuonMinus()->GetProcessManager();
  if(pManager) pManager->RemoveProcess(&fMuMinusIonisation);
  if(pManager) pManager->RemoveProcess(&fMuMinusBremsstrahlung);
  if(pManager) pManager->RemoveProcess(&fMuMinusPairProduction);
  if(pManager) pManager->RemoveProcess(&fMuMinusMultipleScattering);
  if(pManager) pManager->RemoveProcess(&fMuMinusCaptureAtRest);
  pManager = G4TauPlus::TauPlus()->GetProcessManager();
  if(pManager) pManager->RemoveProcess(&fTauPlusIonisation);
  if(pManager) pManager->RemoveProcess(&fTauPlusMultipleScattering);
  pManager = G4TauMinus::TauMinus()->GetProcessManager();
  if(pManager) pManager->RemoveProcess(&fTauMinusIonisation);
  if(pManager) pManager->RemoveProcess(&fTauMinusMultipleScattering);
  }
}

void MuonPhysics::ConstructParticle()
{
  // Mu
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();

  // Tau
  G4TauMinus::TauMinusDefinition();
  G4TauPlus::TauPlusDefinition();
  G4NeutrinoTau::NeutrinoTauDefinition();
  G4AntiNeutrinoTau::AntiNeutrinoTauDefinition();

}

void MuonPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;

  wasActivated = true;
  // Muon Plus Physics
  pManager = G4MuonPlus::MuonPlus()->GetProcessManager();
   // add processes
  pManager->AddProcess(&fMuPlusIonisation, ordInActive,2, 2);

  pManager->AddDiscreteProcess(&fMuPlusBremsstrahlung);

  pManager->AddDiscreteProcess(&fMuPlusPairProduction);

  pManager->AddProcess(&fMuPlusMultipleScattering);
  pManager->SetProcessOrdering(&fMuPlusMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fMuPlusMultipleScattering, idxPostStep,  1);

  // Muon Minus Physics
  pManager = G4MuonMinus::MuonMinus()->GetProcessManager();
   // add processes
  pManager->AddProcess(&fMuMinusIonisation, ordInActive,2, 2);

  pManager->AddDiscreteProcess(&fMuMinusBremsstrahlung);

  pManager->AddDiscreteProcess(&fMuMinusPairProduction);

  pManager->AddProcess(&fMuMinusMultipleScattering);
  pManager->SetProcessOrdering(&fMuMinusMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fMuMinusMultipleScattering, idxPostStep,  1);
  pManager->AddRestProcess(&fMuMinusCaptureAtRest);

  // Tau Plus Physics
  pManager = G4TauPlus::TauPlus()->GetProcessManager();
   // add processes
  pManager->AddProcess(&fTauPlusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&fTauPlusMultipleScattering);
  pManager->SetProcessOrdering(&fTauPlusMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fTauPlusMultipleScattering, idxPostStep,  1);

  // Tau Minus Physics
  pManager = G4TauMinus::TauMinus()->GetProcessManager();
   // add processes
  pManager->AddProcess(&fTauMinusIonisation, ordInActive,2, 2);

  pManager->AddProcess(&fTauMinusMultipleScattering);
  pManager->SetProcessOrdering(&fTauMinusMultipleScattering, idxAlongStep,  1);
  pManager->SetProcessOrdering(&fTauMinusMultipleScattering, idxPostStep,  1);

}



// 2002 by J.P. Wellisch
