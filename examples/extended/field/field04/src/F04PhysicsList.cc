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
//
/// \file field/field04/src/F04PhysicsList.cc
/// \brief Implementation of the F04PhysicsList class
//

#include "F04PhysicsList.hh"
#include "F04PhysicsListMessenger.hh"

#include "G4StepLimiterPhysics.hh"
#include "G4OpticalPhysics.hh"

#include "G4LossTableManager.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

//#include "G4PhysListFactory.hh"
#include "FTFP_BERT.hh"
#include "QGSP_BERT.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "F04StepMax.hh"

#include "G4ProcessTable.hh"

#include "G4PionDecayMakeSpin.hh"
#include "G4DecayWithSpin.hh"

#include "G4DecayTable.hh"
#include "G4MuonDecayChannelWithSpin.hh"
#include "G4MuonRadiativeDecayChannelWithSpin.hh"

#include "G4MuonMinusCapture.hh"
#include "G4MuMinusCapturePrecompound.hh"

#include "G4SystemOfUnits.hh"
#include "G4AutoDelete.hh"

G4ThreadLocal F04StepMax* F04PhysicsList::fStepMaxProcess = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04PhysicsList::F04PhysicsList(G4String physName) : G4VModularPhysicsList()
{
    G4LossTableManager::Instance();

    defaultCutValue  = 1.*mm;

    fMessenger = new F04PhysicsListMessenger(this);

    SetVerboseLevel(1);

//    G4PhysListFactory factory;
    G4VModularPhysicsList* phys = NULL;
    if (physName == "QGSP_BERT") {
       phys = new QGSP_BERT;
    } else {
       phys = new FTFP_BERT;
    }

//    if (factory.IsReferencePhysList(physName))
//       phys =factory.GetReferencePhysList(physName);

    // Physics List is defined via environment variable PHYSLIST
//    if (!phys) phys = factory.ReferencePhysList();

    if (!phys) G4Exception("F04PhysicsList::F04PhysicsList","InvalidSetup",
                              FatalException,"PhysicsList does not exist");

    for (G4int i = 0; ; ++i) {
       G4VPhysicsConstructor* elem =
                  const_cast<G4VPhysicsConstructor*> (phys->GetPhysics(i));
       if (elem == NULL) break;
       G4cout << "RegisterPhysics: " << elem->GetPhysicsName() << G4endl;
       RegisterPhysics(elem);
    }

    RegisterPhysics(new G4StepLimiterPhysics());
    RegisterPhysics(new G4OpticalPhysics());

    fMaxChargedStep = DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04PhysicsList::~F04PhysicsList()
{
    delete fMessenger;

    //delete fStepMaxProcess;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04PhysicsList::ConstructParticle()
{
    G4VModularPhysicsList::ConstructParticle();

    G4GenericIon::GenericIonDefinition();

    G4DecayTable* muonPlusDecayTable = new G4DecayTable();
    muonPlusDecayTable -> Insert(new
                           G4MuonDecayChannelWithSpin("mu+",0.986));
    muonPlusDecayTable -> Insert(new
                           G4MuonRadiativeDecayChannelWithSpin("mu+",0.014));
    G4MuonPlus::MuonPlusDefinition() -> SetDecayTable(muonPlusDecayTable);

    G4DecayTable* muonMinusDecayTable = new G4DecayTable();
    muonMinusDecayTable -> Insert(new
                            G4MuonDecayChannelWithSpin("mu-",0.986));
    muonMinusDecayTable -> Insert(new
                            G4MuonRadiativeDecayChannelWithSpin("mu-",0.014));
    G4MuonMinus::MuonMinusDefinition() -> SetDecayTable(muonMinusDecayTable);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04PhysicsList::ConstructProcess()
{
    G4VModularPhysicsList::ConstructProcess();

    fStepMaxProcess = new F04StepMax();
    G4AutoDelete::Register(fStepMaxProcess);
    
    G4DecayWithSpin* decayWithSpin = new G4DecayWithSpin();

    G4ProcessTable* processTable = G4ProcessTable::GetProcessTable();

    G4VProcess* decay;
    decay = processTable->FindProcess("Decay",G4MuonPlus::MuonPlus());

    G4ProcessManager* pmanager;
    pmanager = G4MuonPlus::MuonPlus()->GetProcessManager();

    if (pmanager) {
      if (decay) pmanager->RemoveProcess(decay);
      pmanager->AddProcess(decayWithSpin);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(decayWithSpin, idxPostStep);
      pmanager ->SetProcessOrdering(decayWithSpin, idxAtRest);
    }

    decay = processTable->FindProcess("Decay",G4MuonMinus::MuonMinus());

    pmanager = G4MuonMinus::MuonMinus()->GetProcessManager();

    if (pmanager) {
      if (decay) pmanager->RemoveProcess(decay);
      pmanager->AddProcess(decayWithSpin);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(decayWithSpin, idxPostStep);
      pmanager ->SetProcessOrdering(decayWithSpin, idxAtRest);
    }

    G4VProcess* process = processTable->
         FindProcess("muMinusCaptureAtRest",G4MuonMinus::MuonMinus());

    if (pmanager) {
       if (process) pmanager->RemoveProcess(process);
       process = new G4MuonMinusCapture(new G4MuMinusCapturePrecompound());
       pmanager->AddRestProcess(process);
    }

    G4PionDecayMakeSpin* poldecay = new G4PionDecayMakeSpin();

    decay = processTable->FindProcess("Decay",G4PionPlus::PionPlus());

    pmanager = G4PionPlus::PionPlus()->GetProcessManager();

    if (pmanager) {
      if (decay) pmanager->RemoveProcess(decay);
      pmanager->AddProcess(poldecay);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(poldecay, idxPostStep);
      pmanager ->SetProcessOrdering(poldecay, idxAtRest);
    }

    decay = processTable->FindProcess("Decay",G4PionMinus::PionMinus());

    pmanager = G4PionMinus::PionMinus()->GetProcessManager();

    if (pmanager) {
      if (decay) pmanager->RemoveProcess(decay);
      pmanager->AddProcess(poldecay);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(poldecay, idxPostStep);
      pmanager ->SetProcessOrdering(poldecay, idxAtRest);
    }

    AddStepMax();
}

/*
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04PhysicsList::RemoveFromPhysicsList(const G4String& name)
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
       std::ostringstream message;
       message << "PhysicsList::RemoveFromPhysicsList "<< name << "not found";
       G4Exception(message.str().c_str());
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04PhysicsList::ClearPhysics()
{
    for (G4PhysConstVector::iterator p  = physicsVector->begin();
                                     p != physicsVector->end(); ++p) {
        delete (*p);
    }
    physicsVector->clear();
}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04PhysicsList::SetStepMax(G4double step)
{
  fMaxChargedStep = step ;
  fStepMaxProcess->SetStepMax(fMaxChargedStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F04StepMax* F04PhysicsList::GetStepMaxProcess()
{
  return fStepMaxProcess;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F04PhysicsList::AddStepMax()
{
  // Step limitation seen as a process

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)()){
      G4ParticleDefinition* particle = particleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();

      if (fStepMaxProcess->IsApplicable(*particle) && !particle->IsShortLived())
      {
         if (pmanager) pmanager ->AddDiscreteProcess(fStepMaxProcess);
      }
  }
}
