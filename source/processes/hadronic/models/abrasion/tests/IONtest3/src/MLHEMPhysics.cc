////////////////////////////////////////////////////////////////////////////////
//
#include "MLHEMPhysics.hh"

#include "globals.hh"
#include "G4ios.hh"
#include <iomanip>
////////////////////////////////////////////////////////////////////////////////
//
MLHEMPhysics::MLHEMPhysics (const G4String& name)
  : G4VPhysicsConstructor(name), mode(name)
{}
////////////////////////////////////////////////////////////////////////////////
//
MLHEMPhysics::~MLHEMPhysics ()
{}
////////////////////////////////////////////////////////////////////////////////
//
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4Proton.hh"
#include "G4AntiProton.hh"


#include "G4PionPlus.hh"
#include "G4PionMinus.hh"

#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"

#include "G4SigmaPlus.hh"
#include "G4SigmaMinus.hh"
#include "G4XiMinus.hh"
#include "G4OmegaMinus.hh"

#include "G4AntiSigmaPlus.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4AntiXiMinus.hh"
#include "G4AntiOmegaMinus.hh"


////////////////////////////////////////////////////////////////////////////////
//
void MLHEMPhysics::ConstructParticle ()
{
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();

  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();

  G4SigmaPlus::SigmaPlusDefinition();
  G4AntiSigmaPlus::AntiSigmaPlusDefinition();
  G4SigmaMinus::SigmaMinusDefinition();
  G4AntiSigmaMinus::AntiSigmaMinusDefinition();
  G4XiMinus::XiMinusDefinition();
  G4AntiXiMinus::AntiXiMinusDefinition();
  G4OmegaMinus::OmegaMinusDefinition();
  G4AntiOmegaMinus::AntiOmegaMinusDefinition();
}
////////////////////////////////////////////////////////////////////////////////

#include "G4ProcessManager.hh"

#include "G4hIonisation.hh"
#include "G4hLowEnergyIonisation.hh"

void MLHEMPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;

  // first to decide what prcesses shall be included according to 
  // the mode string passed on
  //

  G4bool STEM = true ;

  if (mode == "l.e. EM") {
    STEM = false;
  }
  // PionPlus
  pManager = G4PionPlus::PionPlus()->GetProcessManager();
  // add process
  if (STEM) {
    G4hIonisation* thePionPlusIonisation = new G4hIonisation();
    pManager->AddProcess(thePionPlusIonisation, ordInActive,2, 2);

  } else {
    G4hLowEnergyIonisation* thePionPlusIonisation = new G4hLowEnergyIonisation();
    pManager->AddProcess(thePionPlusIonisation, ordInActive,2, 2);
  }

  pManager->AddProcess(&thePionPlusMult);
  pManager->SetProcessOrdering(&thePionPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&thePionPlusMult, idxPostStep, 1);

  // PionMinus
  pManager = G4PionMinus::PionMinus()->GetProcessManager();

  if (STEM) {
    G4hIonisation* thePionMinusIonisation = new G4hIonisation();
    pManager->AddProcess(thePionMinusIonisation, ordInActive,2, 2);
  } else {
    G4hLowEnergyIonisation* thePionMinusIonisation = new G4hLowEnergyIonisation();
    pManager->AddProcess(thePionMinusIonisation, ordInActive,2, 2);
  }
  pManager->AddProcess(&thePionMinusMult);
  pManager->SetProcessOrdering(&thePionMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&thePionMinusMult, idxPostStep, 1);

  // KaonPlus
  pManager = G4KaonPlus::KaonPlus()->GetProcessManager();
  // add process
  if (STEM) {
    G4hIonisation* theKaonPlusIonisation = new G4hIonisation();
    pManager->AddProcess(theKaonPlusIonisation, ordInActive,2, 2);
  } else {
    G4hLowEnergyIonisation* theKaonPlusIonisation = new G4hLowEnergyIonisation();
    pManager->AddProcess(theKaonPlusIonisation, ordInActive,2, 2);
  }
  pManager->AddProcess(&theKaonPlusMult);
  pManager->SetProcessOrdering(&theKaonPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theKaonPlusMult, idxPostStep, 1);

  // KaonMinus
  pManager = G4KaonMinus::KaonMinus()->GetProcessManager();
  // add process
  if (STEM) {
    G4hIonisation* theKaonMinusIonisation = new G4hIonisation();
    pManager->AddProcess(theKaonMinusIonisation, ordInActive,2, 2);
  } else {
    G4hLowEnergyIonisation* theKaonMinusIonisation = new G4hLowEnergyIonisation();
    pManager->AddProcess(theKaonMinusIonisation, ordInActive,2, 2);
  }
  pManager->AddProcess(&theKaonMinusMult);
  pManager->SetProcessOrdering(&theKaonMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theKaonMinusMult, idxPostStep, 1);

  // Proton
  pManager = G4Proton::Proton()->GetProcessManager();
  // add process
  if (STEM) {
    G4hIonisation* theProtonIonisation = new G4hIonisation();
    pManager->AddProcess(theProtonIonisation, ordInActive,2, 2);
  } else {
    G4hLowEnergyIonisation* theProtonIonisation = new G4hLowEnergyIonisation();
    pManager->AddProcess(theProtonIonisation, ordInActive,2, 2);
  }
  pManager->AddProcess(&theProtonMult);
  pManager->SetProcessOrdering(&theProtonMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theProtonMult, idxPostStep, 1);

  // anti-Proton
  pManager = G4AntiProton::AntiProton()->GetProcessManager();
  // add process
  if (STEM) {
    G4hIonisation* theAntiProtonIonisation = new G4hIonisation();
    pManager->AddProcess(theAntiProtonIonisation, ordInActive,2, 2);
  } else {
    G4hLowEnergyIonisation* theAntiProtonIonisation = new G4hLowEnergyIonisation();
    pManager->AddProcess(theAntiProtonIonisation, ordInActive,2, 2);
  }
  pManager->AddProcess(&theAntiProtonMult);
  pManager->SetProcessOrdering(&theAntiProtonMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theAntiProtonMult, idxPostStep, 1);

  // SigmaMinus
  pManager = G4SigmaMinus::SigmaMinus()->GetProcessManager();
  // add process
  if (STEM) {
    G4hIonisation* theSigmaMinusIonisation = new G4hIonisation();
    pManager->AddProcess(theSigmaMinusIonisation, ordInActive,2, 2);
  } else {
    G4hLowEnergyIonisation* theSigmaMinusIonisation = new G4hLowEnergyIonisation();
    pManager->AddProcess(theSigmaMinusIonisation, ordInActive,2, 2);
  }
  pManager->AddProcess(&theSigmaMinusMult);
  pManager->SetProcessOrdering(&theSigmaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theSigmaMinusMult, idxPostStep, 1);

  // anti-SigmaMinus
  pManager = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();
  // add process
  if (STEM) {
    G4hIonisation* theAntiSigmaMinusIonisation = new G4hIonisation();
    pManager->AddProcess(theAntiSigmaMinusIonisation, ordInActive,2, 2);
  } else {
    G4hLowEnergyIonisation* theAntiSigmaMinusIonisation = new G4hLowEnergyIonisation();
    pManager->AddProcess(theAntiSigmaMinusIonisation, ordInActive,2, 2);
  }
  pManager->AddProcess(&theAntiSigmaMinusMult);
  pManager->SetProcessOrdering(&theAntiSigmaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theAntiSigmaMinusMult, idxPostStep, 1);

  // SigmaPlus
  pManager = G4SigmaPlus::SigmaPlus()->GetProcessManager();
  // add process
  if (STEM) {
    G4hIonisation* theSigmaPlusIonisation = new G4hIonisation();
    pManager->AddProcess(theSigmaPlusIonisation, ordInActive,2, 2);
  } else {
    G4hLowEnergyIonisation*  theSigmaPlusIonisation = new G4hLowEnergyIonisation();
    pManager->AddProcess(theSigmaPlusIonisation, ordInActive,2, 2);
  }  
  pManager->AddProcess(&theSigmaPlusMult);
  pManager->SetProcessOrdering(&theSigmaPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theSigmaPlusMult, idxPostStep, 1);

  // anti-SigmaPlus
  pManager = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();
  // add process
  if (STEM) {
    G4hIonisation* theAntiSigmaPlusIonisation = new G4hIonisation();
    pManager->AddProcess(theAntiSigmaPlusIonisation, ordInActive,2, 2);
  } else {
    G4hLowEnergyIonisation* theAntiSigmaPlusIonisation = new G4hLowEnergyIonisation();
    pManager->AddProcess(theAntiSigmaPlusIonisation, ordInActive,2, 2);
  }
  pManager->AddProcess(&theAntiSigmaPlusMult);
  pManager->SetProcessOrdering(&theAntiSigmaPlusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theAntiSigmaPlusMult, idxPostStep, 1);

  // XiMinus
  pManager = G4XiMinus::XiMinus()->GetProcessManager();
  // add process
  if (STEM) {
    G4hIonisation* theXiMinusIonisation = new G4hIonisation();
    pManager->AddProcess(theXiMinusIonisation, ordInActive,2, 2);
  } else {
    G4hLowEnergyIonisation* theXiMinusIonisation = new G4hLowEnergyIonisation();
    pManager->AddProcess(theXiMinusIonisation, ordInActive,2, 2);
  }
  pManager->AddProcess(&theXiMinusMult);
  pManager->SetProcessOrdering(&theXiMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theXiMinusMult, idxPostStep, 1);

  // anti-XiMinus
  pManager = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();
  // add process
  if (STEM) {
    G4hIonisation* theAntiXiMinusIonisation = new G4hIonisation();
    pManager->AddProcess(theAntiXiMinusIonisation, ordInActive,2, 2);
  } else {
    G4hLowEnergyIonisation* theAntiXiMinusIonisation = new G4hLowEnergyIonisation();
    pManager->AddProcess(theAntiXiMinusIonisation, ordInActive,2, 2);
  }
  pManager->AddProcess(&theAntiXiMinusMult);
  pManager->SetProcessOrdering(&theAntiXiMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theAntiXiMinusMult, idxPostStep, 1);

  // OmegaMinus
  pManager = G4OmegaMinus::OmegaMinus()->GetProcessManager();
  // add process
  if (STEM) {
    G4hIonisation* theOmegaMinusIonisation = new G4hIonisation();
    pManager->AddProcess(theOmegaMinusIonisation, ordInActive,2, 2);
  } else {
    G4hLowEnergyIonisation* theOmegaMinusIonisation = new G4hLowEnergyIonisation();
    pManager->AddProcess(theOmegaMinusIonisation, ordInActive,2, 2);
  }
  pManager->AddProcess(&theOmegaMinusMult);
  pManager->SetProcessOrdering(&theOmegaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theOmegaMinusMult, idxPostStep, 1);

  // anti-OmegaMinus
  pManager = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();
  // add process
  if (STEM) {
    G4hIonisation* theAntiOmegaMinusIonisation = new G4hIonisation();
    pManager->AddProcess(theAntiOmegaMinusIonisation, ordInActive,2, 2);
  } else {
    G4hLowEnergyIonisation* theAntiOmegaMinusIonisation = new G4hLowEnergyIonisation();
    pManager->AddProcess(theAntiOmegaMinusIonisation, ordInActive,2, 2);
  }
  pManager->AddProcess(&theAntiOmegaMinusMult);
  pManager->SetProcessOrdering(&theAntiOmegaMinusMult, idxAlongStep, 1);
  pManager->SetProcessOrdering(&theAntiOmegaMinusMult, idxPostStep, 1);

}
////////////////////////////////////////////////////////////////////////////////
