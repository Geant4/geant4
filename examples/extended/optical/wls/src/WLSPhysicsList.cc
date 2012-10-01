//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file optical/wls/src/WLSPhysicsList.cc
/// \brief Implementation of the WLSPhysicsList class
//
//
//

#include "WLSPhysicsList.hh"
#include "WLSPhysicsListMessenger.hh"

#include "WLSExtraPhysics.hh"
#include "WLSOpticalPhysics.hh"

#include "G4LossTableManager.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4PhysListFactory.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "WLSStepMax.hh"

#include "G4ProcessTable.hh"

#include "G4PionDecayMakeSpin.hh"
#include "G4DecayWithSpin.hh"

#include "G4DecayTable.hh"
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"

#include "G4RadioactiveDecayPhysics.hh"

#include "G4SystemOfUnits.hh"

WLSPhysicsList::WLSPhysicsList(G4String physName) : G4VModularPhysicsList()
{
    G4LossTableManager::Instance();

    defaultCutValue  = 1.*mm;
    fCutForGamma     = defaultCutValue;
    fCutForElectron  = defaultCutValue;
    fCutForPositron  = defaultCutValue;

    G4PhysListFactory factory;
    G4VModularPhysicsList* phys = NULL;
    if (factory.IsReferencePhysList(physName)) {
       phys = factory.GetReferencePhysList(physName);
       if(!phys)G4Exception("WLSPhysicsList::WLSPhysicsList","InvalidSetup",
                            FatalException,"PhysicsList does not exist");
       fMessenger = new WLSPhysicsListMessenger(this);
    }

    for (G4int i = 0; ; ++i) {
       G4VPhysicsConstructor* elem =
                  const_cast<G4VPhysicsConstructor*> (phys->GetPhysics(i));
       if (elem == NULL) break;
       G4cout << "RegisterPhysics: " << elem->GetPhysicsName() << G4endl;
       RegisterPhysics(elem);
    }

    AbsorptionOn = true;

    physicsVector->push_back(new WLSExtraPhysics());
    physicsVector->push_back(opticalPhysics = new WLSOpticalPhysics(AbsorptionOn));

    physicsVector->push_back(new G4RadioactiveDecayPhysics());

    stepMaxProcess = new WLSStepMax();
}

WLSPhysicsList::~WLSPhysicsList()
{
    delete fMessenger;

    delete stepMaxProcess;
}

void WLSPhysicsList::ClearPhysics()
{
    for (G4PhysConstVector::iterator p  = physicsVector->begin();
                                     p != physicsVector->end(); ++p) {
        delete (*p);
    }
    physicsVector->clear();
}

void WLSPhysicsList::ConstructParticle()
{
    G4VModularPhysicsList::ConstructParticle();

    G4DecayTable* MuonPlusDecayTable = new G4DecayTable();
    MuonPlusDecayTable -> Insert(new
                           G4MuonDecayChannelWithSpin("mu+",0.986));
    MuonPlusDecayTable -> Insert(new
                           G4MuonRadiativeDecayChannelWithSpin("mu+",0.014));
    G4MuonPlus::MuonPlusDefinition() -> SetDecayTable(MuonPlusDecayTable);

    G4DecayTable* MuonMinusDecayTable = new G4DecayTable();
    MuonMinusDecayTable -> Insert(new
                            G4MuonDecayChannelWithSpin("mu-",0.986));
    MuonMinusDecayTable -> Insert(new
                            G4MuonRadiativeDecayChannelWithSpin("mu-",0.014));
    G4MuonMinus::MuonMinusDefinition() -> SetDecayTable(MuonMinusDecayTable);

}

void WLSPhysicsList::ConstructProcess()
{
    G4VModularPhysicsList::ConstructProcess();

    SetVerbose(0);

    G4DecayWithSpin* decayWithSpin = new G4DecayWithSpin();

    G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();

    G4VProcess* decay;
    decay = processTable->FindProcess("Decay",G4MuonPlus::MuonPlus());

    G4ProcessManager* fManager;
    fManager = G4MuonPlus::MuonPlus()->GetProcessManager();

    if (fManager) {
      if (decay) fManager->RemoveProcess(decay);
      fManager->AddProcess(decayWithSpin);
      // set ordering for PostStepDoIt and AtRestDoIt
      fManager ->SetProcessOrdering(decayWithSpin, idxPostStep);
      fManager ->SetProcessOrdering(decayWithSpin, idxAtRest);
    }

    decay = processTable->FindProcess("Decay",G4MuonMinus::MuonMinus());

    fManager = G4MuonMinus::MuonMinus()->GetProcessManager();

    if (fManager) {
      if (decay) fManager->RemoveProcess(decay);
      fManager->AddProcess(decayWithSpin);
      // set ordering for PostStepDoIt and AtRestDoIt
      fManager ->SetProcessOrdering(decayWithSpin, idxPostStep);
      fManager ->SetProcessOrdering(decayWithSpin, idxAtRest);
    }

    G4PionDecayMakeSpin* poldecay = new G4PionDecayMakeSpin();

    decay = processTable->FindProcess("Decay",G4PionPlus::PionPlus());

    fManager = G4PionPlus::PionPlus()->GetProcessManager();

    if (fManager) {
      if (decay) fManager->RemoveProcess(decay);
      fManager->AddProcess(poldecay);
      // set ordering for PostStepDoIt and AtRestDoIt
      fManager ->SetProcessOrdering(poldecay, idxPostStep);
      fManager ->SetProcessOrdering(poldecay, idxAtRest);
    }

    decay = processTable->FindProcess("Decay",G4PionMinus::PionMinus());

    fManager = G4PionMinus::PionMinus()->GetProcessManager();

    if (fManager) {
      if (decay) fManager->RemoveProcess(decay);
      fManager->AddProcess(poldecay);
      // set ordering for PostStepDoIt and AtRestDoIt
      fManager ->SetProcessOrdering(poldecay, idxPostStep);
      fManager ->SetProcessOrdering(poldecay, idxAtRest);
    }

    AddStepMax();

}

void WLSPhysicsList::RemoveFromPhysicsList(const G4String& name)
{
    G4bool success = false;
    for (G4PhysConstVector::iterator p  = physicsVector->begin();
                                     p != physicsVector->end(); ++p) {
        G4VPhysicsConstructor* e = (*p);
        if (e->GetPhysicsName() == name) {
           physicsVector->erase(p);
           success = true;
           break;
        }
    }
    if (!success) {
       G4ExceptionDescription message;
       message << "PhysicsList::RemoveFromEMPhysicsList "<< name << "not found";
       G4Exception("example WLSPhysicsList::RemoveFromPhysicsList()",
       "ExamWLSPhysicsList01",FatalException,message);
    }
}

void WLSPhysicsList::SetAbsorption(G4bool toggle)
{
       AbsorptionOn = toggle;
       RemoveFromPhysicsList("Optical");
       physicsVector->push_back(opticalPhysics = new WLSOpticalPhysics(toggle));
       opticalPhysics->ConstructProcess();
}

void WLSPhysicsList::SetCuts()
{
    if (verboseLevel >0) {
        G4cout << "WLSPhysicsList::SetCuts:";
        G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length")
               << G4endl;
    }

    // set cut values for gamma at first and for e- second and next for e+,
    // because some processes for e+/e- need cut values for gamma
    SetCutValue(fCutForGamma, "gamma");
    SetCutValue(fCutForElectron, "e-");
    SetCutValue(fCutForPositron, "e+");

    if (verboseLevel>0) DumpCutValuesTable();
}

void WLSPhysicsList::SetCutForGamma(G4double cut)
{
    fCutForGamma = cut;
    SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

void WLSPhysicsList::SetCutForElectron(G4double cut)
{
    fCutForElectron = cut;
    SetParticleCuts(fCutForElectron, G4Electron::Electron());
}

void WLSPhysicsList::SetCutForPositron(G4double cut)
{
    fCutForPositron = cut;
    SetParticleCuts(fCutForPositron, G4Positron::Positron());
}

void WLSPhysicsList::SetStepMax(G4double step)
{
  MaxChargedStep = step;
  stepMaxProcess->SetStepMax(MaxChargedStep);
}

WLSStepMax* WLSPhysicsList::GetStepMaxProcess()
{
  return stepMaxProcess;
}

void WLSPhysicsList::AddStepMax()
{
  // Step limitation seen as a process

  theParticleIterator->reset();
  while ((*theParticleIterator)()){
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();

      if (stepMaxProcess->IsApplicable(*particle) && !particle->IsShortLived())
      {
         if (pmanager) pmanager ->AddDiscreteProcess(stepMaxProcess);
      }
  }
}

void WLSPhysicsList::SetNbOfPhotonsCerenkov(G4int MaxNumber)
{
   opticalPhysics->SetNbOfPhotonsCerenkov(MaxNumber);
}

void WLSPhysicsList::SetVerbose(G4int verbose)
{
   opticalPhysics->GetCerenkovProcess()->SetVerboseLevel(verbose);
   opticalPhysics->GetScintillationProcess()->SetVerboseLevel(verbose);
   opticalPhysics->GetAbsorptionProcess()->SetVerboseLevel(verbose);
   opticalPhysics->GetRayleighScatteringProcess()->SetVerboseLevel(verbose);
   opticalPhysics->GetMieHGScatteringProcess()->SetVerboseLevel(verbose);
   opticalPhysics->GetBoundaryProcess()->SetVerboseLevel(verbose);
}
