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
// $Id: PhysicsList.cc,v 1.1 2007-05-24 21:57:03 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, May 2007. Copy of PhysicsList found in 
//                         examples/extended/electromagnetic/TestEm7. 
//                         Also implemented dumpInfo.
//
#include "PhysicsList.hh"

#include "G4LossTableManager.hh"
#include "G4UnitsTable.hh"
#include "LowEnergyPhysicsList.hh" 
#include "PhysicsListMessenger.hh"
#include "StandardPhysicsList.hh"

// Bosons
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"

// leptons
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

// Mesons
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionZero.hh"
#include "G4Eta.hh"
#include "G4EtaPrime.hh"

#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonZero.hh"
#include "G4AntiKaonZero.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"

// Baryons
#include "G4Proton.hh"
#include "G4AntiProton.hh"
#include "G4Neutron.hh"
#include "G4AntiNeutron.hh"

// Nuclei
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"

#include "G4Decay.hh"
#include "G4ProcessManager.hh"
#include "G4ProductionCutsTable.hh"
#include "ConfigData.hh"

#include "G4EmProcessOptions.hh"

PhysicsList::PhysicsList() 
{
  G4LossTableManager::Instance();
  
  currentDefaultCut   = 1.0*mm;
  cutForGamma         = currentDefaultCut;
  cutForElectron      = currentDefaultCut;
  cutForPositron      = currentDefaultCut;

  pMessenger = new PhysicsListMessenger(this);

  SetVerboseLevel(1);
}

PhysicsList::~PhysicsList()
{
  delete pMessenger;
}

void PhysicsList::ConstructParticle()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  
  // gamma
  G4Gamma::GammaDefinition();
  
  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
  
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();  
  
  // mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();

  // barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();

  // ions
  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4Alpha::AlphaDefinition();
  G4GenericIon::GenericIonDefinition();
}


void PhysicsList::ConstructProcess()
{
  // Transportation
  AddTransportation();
  
  pImpl->ConstructProcess();
  
  G4Decay* fDecayProcess = new G4Decay();
  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    
    if (fDecayProcess->IsApplicable(*particle) && !particle->IsShortLived()) { 
      
      pmanager ->AddProcess(fDecayProcess);
      
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(fDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(fDecayProcess, idxAtRest);
      
    }
  }

  DumpInfo();
}

void PhysicsList::AddPhysicsList(const G4String& name)
{
  if (verboseLevel>0) {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }
  if (name == "Standard") {
    pImpl = new StandardPhysicsList();
      
  } else if (name == "LowEnergy") {
    pImpl = new LowEnergyPhysicsList();
  } 
  else {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
           << " is not defined"
           << G4endl;
    assert(0);
  }
}


void PhysicsList::SetCuts()
{    
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
}

void PhysicsList::SetCutForGamma(G4double cut)
{
  cutForGamma = cut;
  SetParticleCuts(cutForGamma, G4Gamma::Gamma());
}

void PhysicsList::SetCutForElectron(G4double cut)
{
  cutForElectron = cut;
  SetParticleCuts(cutForElectron, G4Electron::Electron());
}

void PhysicsList::SetCutForPositron(G4double cut)
{
  cutForPositron = cut;
  SetParticleCuts(cutForPositron, G4Positron::Positron());
}

void PhysicsList::DumpInfo()
{
  std::ofstream& log = ConfigData::GetConfigFile();
  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    log << "===================================================================" << G4endl;
    log << "Processes for particle "<<particle->GetParticleName()<<G4endl;
    
    
    // loop over all processes
    G4ProcessVector* theProcessList = pmanager->GetProcessList();
    
    for (G4int idx=0; idx <theProcessList->entries(); idx++){
      // process name/type
      log << "[" << idx << "]";
      log << "=== process[" << ((*theProcessList)(idx))->GetProcessName()<< " :";
      log << G4VProcess::GetProcessTypeName( ((*theProcessList)(idx))->GetProcessType() )<< "]"<<G4endl;      
      log << "  Ordering::     ";
      log << "        AtRest             AlongStep          PostStep   ";
      log << G4endl;
      log << "                 ";
      log << "   GetPIL/    DoIt    GetPIL/    DoIt    GetPIL/    DoIt ";
      log << G4endl;   
      log << "  Ordering::      " << G4endl;
      G4VProcess* proc = (*theProcessList)(idx);
      log << "  index           ";
      log<<std::setw(9)<<pmanager->GetAtRestIndex(proc, typeGPIL)<<std::setw(9)<<pmanager->GetAtRestIndex(proc, typeDoIt)<<std::setw(9)<<
	pmanager->GetAlongStepIndex(proc, typeGPIL)<<std::setw(9)<<pmanager->GetAlongStepIndex(proc, typeDoIt)<<std::setw(9)<<
	pmanager->GetPostStepIndex(proc, typeGPIL)<<std::setw(9)<<pmanager->GetPostStepIndex(proc, typeDoIt)<<G4endl;
      
    }
  }
}
