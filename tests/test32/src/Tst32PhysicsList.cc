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
#include "globals.hh"
#include "Tst32PhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4ios.hh"
#include <iomanip>                

#include "G4Mars5GeV.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst32PhysicsList::Tst32PhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 1.0*cm;
   SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst32PhysicsList::~Tst32PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst32PhysicsList::ConstructParticle()
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst32PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst32PhysicsList::ConstructLeptons()
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst32PhysicsList::ConstructMesons()
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst32PhysicsList::ConstructBaryons()
{
  //  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();

  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
  
  

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst32PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructMarsProc();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

////#include "G4UserSpecialCuts.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst32PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
    //electron
      G4VProcess* theeminusMultipleScattering = new G4eMultipleScattering();
      G4VProcess* theeminusIonisation         = new G4eIonisation();
      G4VProcess* theeminusBremsstrahlung     = new G4eBremsstrahlung();
      //
      // add processes
      pmanager->AddProcess(theeminusMultipleScattering);
      pmanager->AddProcess(theeminusIonisation);
      pmanager->AddProcess(theeminusBremsstrahlung);
      //      
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(theeminusIonisation,         idxAlongStep,2);
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(theeminusIonisation,         idxPostStep,2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung,     idxPostStep,3);

    } else if (particleName == "e+") {
    //positron
      G4VProcess* theeplusMultipleScattering = new G4eMultipleScattering();
      G4VProcess* theeplusIonisation         = new G4eIonisation();
      G4VProcess* theeplusBremsstrahlung     = new G4eBremsstrahlung();
      G4VProcess* theeplusAnnihilation       = new G4eplusAnnihilation();
      //
      // add processes
      pmanager->AddProcess(theeplusMultipleScattering);
      pmanager->AddProcess(theeplusIonisation);
      pmanager->AddProcess(theeplusBremsstrahlung);
      pmanager->AddProcess(theeplusAnnihilation);
      //
      // set ordering for AtRestDoIt
      pmanager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);
      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(theeplusIonisation,         idxAlongStep,2);
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(theeplusIonisation,         idxPostStep,2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung,     idxPostStep,3);
      pmanager->SetProcessOrdering(theeplusAnnihilation,       idxPostStep,4);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
    //muon  
      G4VProcess* aMultipleScattering = new G4MuMultipleScattering();
      G4VProcess* aBremsstrahlung     = new G4MuBremsstrahlung();
      G4VProcess* aPairProduction     = new G4MuPairProduction();
      G4VProcess* anIonisation        = new G4MuIonisation();
      //
      // add processes
      pmanager->AddProcess(anIonisation);
      pmanager->AddProcess(aMultipleScattering);
      pmanager->AddProcess(aBremsstrahlung);
      pmanager->AddProcess(aPairProduction);
      //
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(anIonisation,        idxAlongStep,2);
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(anIonisation,        idxPostStep,2);
      pmanager->SetProcessOrdering(aBremsstrahlung,     idxPostStep,3);
      pmanager->SetProcessOrdering(aPairProduction,     idxPostStep,4);

    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
     // all others charged particles except geantino
     G4VProcess* aMultipleScattering = new G4hMultipleScattering();
     G4VProcess* anIonisation        = new G4hIonisation();     
     ////G4VProcess*  theUserCuts = new G4UserSpecialCuts();
     
     //
     // add processes
     pmanager->AddProcess(anIonisation);
     pmanager->AddProcess(aMultipleScattering);    
     ////pmanager->AddProcess(theUserCuts);
     
     //
     // set ordering for AlongStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,1);
     pmanager->SetProcessOrdering(anIonisation,        idxAlongStep,2);
     
     //
     // set ordering for PostStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep,1);
     pmanager->SetProcessOrdering(anIonisation,        idxPostStep,2);
     ////pmanager->SetProcessOrdering(theUserCuts,     idxPostStep,3);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4Decay.hh"
void Tst32PhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
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
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
#include "Tst32EminCut.hh"

void Tst32PhysicsList::ConstructMarsProc()
{
  Tst32EminCut* killProcess = new Tst32EminCut();
  killProcess->SetEminThreshold(1.0 *MeV); 

  // for neutron/proton/anti_proton
  theParticleIterator->reset();

  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager * pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    // baryon
    if (particleName == "neutron") {
      G4HadronInelasticProcess * theParticleIEProc = 
	new G4NeutronInelasticProcess();
      G4Mars5GeV * theParticleIE = new G4Mars5GeV();
      theParticleIE->SetMinEnergy(0.*MeV);
      theParticleIEProc->RegisterMe( theParticleIE );
      pmanager->AddDiscreteProcess( theParticleIEProc );

      pmanager->AddDiscreteProcess( killProcess);
    } 
    else if (particleName == "proton") {
      G4HadronInelasticProcess * theParticleIEProc = 
	new G4ProtonInelasticProcess();
      G4Mars5GeV * theParticleIE = new G4Mars5GeV();
      theParticleIE->SetMinEnergy(0.*MeV);
      theParticleIEProc->RegisterMe( theParticleIE );
      pmanager->AddDiscreteProcess( theParticleIEProc );

      pmanager->AddDiscreteProcess( killProcess);
    }
    else if (particleName == "antiproton") {
      G4HadronInelasticProcess * theParticleIEProc = 
	new G4AntiProtonInelasticProcess();
      G4Mars5GeV * theParticleIE = new G4Mars5GeV();
      theParticleIE->SetMinEnergy(0.*MeV);
      theParticleIEProc->RegisterMe( theParticleIE );
      pmanager->AddDiscreteProcess( theParticleIEProc );

      pmanager->AddDiscreteProcess( killProcess);
    }
    // meson
    else if (particleName == "pi-") {
      G4HadronInelasticProcess * theParticleIEProc = 
	new G4PionMinusInelasticProcess();
      G4Mars5GeV * theParticleIE = new G4Mars5GeV();
      theParticleIE->SetMinEnergy(0.*MeV);
      theParticleIEProc->RegisterMe( theParticleIE );
      pmanager->AddDiscreteProcess( theParticleIEProc );

      pmanager->AddDiscreteProcess( killProcess);
    }      
    else if (particleName == "pi+") {
      G4HadronInelasticProcess * theParticleIEProc = 
	new G4PionPlusInelasticProcess();
      G4Mars5GeV * theParticleIE = new G4Mars5GeV();
      theParticleIE->SetMinEnergy(0.*MeV);
      theParticleIEProc->RegisterMe( theParticleIE );
      pmanager->AddDiscreteProcess( theParticleIEProc );

      pmanager->AddDiscreteProcess( killProcess);
    }      
    else if (particleName == "kaon-") {
      G4HadronInelasticProcess * theParticleIEProc = 
	new G4KaonMinusInelasticProcess();
      G4Mars5GeV * theParticleIE = new G4Mars5GeV();
      theParticleIE->SetMinEnergy(0.*MeV);
      theParticleIEProc->RegisterMe( theParticleIE );
      pmanager->AddDiscreteProcess( theParticleIEProc );

      pmanager->AddDiscreteProcess( killProcess);
    }      
    else if (particleName == "kaon+") {
      G4HadronInelasticProcess * theParticleIEProc = 
	new G4KaonPlusInelasticProcess();
      G4Mars5GeV * theParticleIE = new G4Mars5GeV();
      theParticleIE->SetMinEnergy(0.*MeV);
      theParticleIEProc->RegisterMe( theParticleIE );
      pmanager->AddDiscreteProcess( theParticleIEProc );

      pmanager->AddDiscreteProcess( killProcess);
    } 
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst32PhysicsList::SetCuts()
{
  //G4VUserPhysicsList::SetCutsWithDefault method sets 
  //the default cut value for all particle types 
  //
  SetCutsWithDefault();
 
  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....






