// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN02PhysicsList.cc,v 1.4 1999-12-15 14:49:21 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This is a version for maximum particle set
// ------------------------------------------------------------

#include "globals.hh"
#include "ExN02PhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include "g4std/iomanip"                


ExN02PhysicsList::ExN02PhysicsList():  G4VUserPhysicsList(),
 thePhotoElectricEffect(NULL),theComptonScattering(NULL),
 theGammaConversion(NULL),
 theeminusMultipleScattering(NULL),theeminusIonisation(NULL),
 theeminusBremsstrahlung(NULL),
 theeplusMultipleScattering(NULL),theeplusIonisation(NULL),
 theeplusBremsstrahlung(NULL),
 theeplusAnnihilation(NULL)
{
  defaultCutValue = 2.0*mm;
   SetVerboseLevel(1);
}

ExN02PhysicsList::~ExN02PhysicsList()
{
}

void ExN02PhysicsList::ConstructParticle()
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

void ExN02PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();
}

void ExN02PhysicsList::ConstructLeptons()
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

void ExN02PhysicsList::ConstructMesons()
{
  //  mesons
  //    light mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4RhoZero::RhoZeroDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
}

void ExN02PhysicsList::ConstructBaryons()
{
  //  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();

  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}


void ExN02PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

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
///#include "G4UserSpecialCuts.hh"

void ExN02PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
      thePhotoElectricEffect = new G4PhotoElectricEffect();
      theComptonScattering = new G4ComptonScattering();
      theGammaConversion = new G4GammaConversion();
      // add processes  
      pmanager->AddDiscreteProcess(thePhotoElectricEffect);
      pmanager->AddDiscreteProcess(theComptonScattering);
      pmanager->AddDiscreteProcess(theGammaConversion);
      
    } else if (particleName == "e-") {
    //electron
      // Construct processes for electron
      theeminusMultipleScattering = new G4MultipleScattering();
      theeminusIonisation = new G4eIonisation();
      theeminusBremsstrahlung = new G4eBremsstrahlung();
      // add processes
      pmanager->AddProcess(theeminusMultipleScattering);
      pmanager->AddProcess(theeminusIonisation);
      pmanager->AddProcess(theeminusBremsstrahlung);      
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxAlongStep,  1);
      pmanager->SetProcessOrdering(theeminusIonisation, idxAlongStep,  2);
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(theeminusIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung, idxPostStep, 3);

    } else if (particleName == "e+") {
    //positron
      // Construct processes for positron
      theeplusMultipleScattering = new G4MultipleScattering();
      theeplusIonisation = new G4eIonisation();
      theeplusBremsstrahlung = new G4eBremsstrahlung();
      theeplusAnnihilation = new G4eplusAnnihilation();
      // add processes
      pmanager->AddProcess(theeplusMultipleScattering);
      pmanager->AddProcess(theeplusIonisation);
      pmanager->AddProcess(theeplusBremsstrahlung);
      pmanager->AddProcess(theeplusAnnihilation);
      // set ordering for AtRestDoIt
      pmanager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxAlongStep,  1);
      pmanager->SetProcessOrdering(theeplusIonisation, idxAlongStep,  2);
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(theeplusIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung, idxPostStep, 3);
      pmanager->SetProcessOrdering(theeplusAnnihilation, idxPostStep, 4);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
    //muon  
     // Construct processes for muon+
     G4VProcess* aMultipleScattering = new G4MultipleScattering();
     G4VProcess* aBremsstrahlung = new G4MuBremsstrahlung();
     G4VProcess* aPairProduction = new G4MuPairProduction();
     G4VProcess* anIonisation = new G4MuIonisation();
      // add processes
     pmanager->AddProcess(anIonisation);
     pmanager->AddProcess(aMultipleScattering);
     pmanager->AddProcess(aBremsstrahlung);
     pmanager->AddProcess(aPairProduction);
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,  1);
      pmanager->SetProcessOrdering(anIonisation, idxAlongStep,  2);
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(anIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(aBremsstrahlung, idxPostStep, 3);
      pmanager->SetProcessOrdering(aPairProduction, idxPostStep, 4);

    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
     // all others charged particles except geantino
     G4VProcess* aMultipleScattering = new G4MultipleScattering();
     G4VProcess* anIonisation = new G4hIonisation();
     //>>G4VProcess*  theUserSpecialCuts = new G4UserSpecialCuts();
     // add processes
     pmanager->AddProcess(anIonisation);
     pmanager->AddProcess(aMultipleScattering);
     /// pmanager->AddProcess(theUserSpecialCuts);
     // set ordering for AtRestDoIt
     //>>pmanager->SetProcessOrderingToFirst(theUserSpecialCuts, idxAtRest);
     // set ordering for AlongStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,  1);
     pmanager->SetProcessOrdering(anIonisation, idxAlongStep,  2);
     // set ordering for PostStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);
     pmanager->SetProcessOrdering(anIonisation, idxPostStep, 2);
     //>>pmanager->SetProcessOrdering(theUserSpecialCuts, idxPostStep, 3);
    }
  }
}

#include "G4Decay.hh"
void ExN02PhysicsList::ConstructGeneral()
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

void ExN02PhysicsList::SetCuts()
{
  if (verboseLevel >1){
    G4cout << "ExN02PhysicsList::SetCuts:";
  }  
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   

}
