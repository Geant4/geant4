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
// $Id: ExTGPhysicsListWithParallel.cc,v 1.1 2009-05-15 16:39:14 arce Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ---------------------------------------------------------------------------

#include "globals.hh"
#include "ExTGPhysicsListWithParallel.hh"

#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"

// ---------------------------------------------------------------------------

ExTGPhysicsListWithParallel::ExTGPhysicsListWithParallel():  G4VUserPhysicsList()
{
  defaultCutValue = 1.0*cm;
   SetVerboseLevel(1);
}

// ---------------------------------------------------------------------------

ExTGPhysicsListWithParallel::~ExTGPhysicsListWithParallel()
{
}

// ---------------------------------------------------------------------------

void ExTGPhysicsListWithParallel::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
}

// ---------------------------------------------------------------------------

void ExTGPhysicsListWithParallel::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();
}

// ---------------------------------------------------------------------------

void ExTGPhysicsListWithParallel::ConstructLeptons()
{
  // leptons
  //  e+/-
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  // mu+/-
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  // nu_e
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  // nu_mu
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

// ---------------------------------------------------------------------------

void ExTGPhysicsListWithParallel::ConstructMesons()
{
  //  mesons
  //    light mesons
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
}

// ---------------------------------------------------------------------------

void ExTGPhysicsListWithParallel::ConstructBaryons()
{
  //  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();

  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}

// ---------------------------------------------------------------------------

void ExTGPhysicsListWithParallel::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

// ---------------------------------------------------------------------------

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

#include "G4StepLimiter.hh"
#include "G4UserSpecialCuts.hh"

// ---------------------------------------------------------------------------

void ExTGPhysicsListWithParallel::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
      // gamma         
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);
      
    } else if (particleName == "e-") {
      //electron
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4eIonisation,       -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,   -1, 3,3);      

    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4eIonisation,       -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,   -1, 3,3);
      pmanager->AddProcess(new G4eplusAnnihilation,  0,-1,4);

    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon  
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4MuIonisation,      -1, 2,2);
      pmanager->AddProcess(new G4MuBremsstrahlung,  -1, 3,3);
      pmanager->AddProcess(new G4MuPairProduction,  -1, 4,4);       
     
    } else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) && 
               (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4hIonisation,       -1, 2,2);
      //step limit
      pmanager->AddProcess(new G4StepLimiter,       -1,-1,3);         
      ///pmanager->AddProcess(new G4UserSpecialCuts,   -1,-1,4);  
    }
  }
}

// ---------------------------------------------------------------------------

#include "G4Decay.hh"
#include "G4ParallelWorldScoringProcess.hh"
#include "G4tgbParallelGeomMgr.hh"
#include "G4VUserParallelWorld.hh"
#include "G4UIcommand.hh"

void ExTGPhysicsListWithParallel::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  G4ParallelWorldScoringProcess* theParallelWorldScoringProcess
    = new G4ParallelWorldScoringProcess("ParaWorldScoringProc");
  theParallelWorldScoringProcess->SetParallelWorld("ParallelScoringWorld");

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
    // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
    pmanager->AddProcess(theParallelWorldScoringProcess);
    pmanager->SetProcessOrderingToLast(theParallelWorldScoringProcess, idxAtRest);
    pmanager->SetProcessOrdering(theParallelWorldScoringProcess, idxAlongStep, 1);
    pmanager->SetProcessOrderingToLast(theParallelWorldScoringProcess, idxPostStep);

  }

  std::vector<G4VUserParallelWorld*> parallelWorlds = G4tgbParallelGeomMgr::GetInstance()->GetParallelWorlds();
  for( size_t ii = 0; ii < parallelWorlds.size(); ii++ ) {
    G4String pwName = parallelWorlds[ii]->GetName();
    G4int index = G4UIcommand::ConvertToInt( pwName.substr(
		   pwName.rfind("_", pwName.length())+1, pwName.length() ).c_str() );
    G4cout << " RegisterParallelWorld " << pwName << " index " << index << " " << G4endl;
    G4tgbParallelGeomMgr::GetInstance()->BuildPhysicsProcess(parallelWorlds[ii]->GetName(), index );
    
  }

}

// ---------------------------------------------------------------------------

void ExTGPhysicsListWithParallel::SetCuts()
{
  //G4VUserPhysicsList::SetCutsWithDefault method sets 
  //the default cut value for all particle types 
  //
  SetCutsWithDefault();
     
  if (verboseLevel>0) DumpCutValuesTable();
}

// ---------------------------------------------------------------------------

