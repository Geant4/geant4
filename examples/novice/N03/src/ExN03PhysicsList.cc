// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: ExN03PhysicsList.cc,v 1.8 2000-05-24 09:53:05 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "ExN03PhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4ios.hh"              

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExN03PhysicsList::ExN03PhysicsList():  G4VUserPhysicsList()
{
  currentDefaultCut = defaultCutValue = 1.0*mm;
  cutForGamma       = defaultCutValue;
  cutForElectron    = defaultCutValue;
  cutForProton      = defaultCutValue;

 SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

ExN03PhysicsList::~ExN03PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExN03PhysicsList::ConstructParticle()
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

void ExN03PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExN03PhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExN03PhysicsList::ConstructMesons()
{
 //  mesons
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

void ExN03PhysicsList::ConstructBaryons()
{
//  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExN03PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void ExN03PhysicsList::ConstructEM()
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
      G4VProcess* theeminusMultipleScattering = new G4MultipleScattering();
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
      G4VProcess* theeplusMultipleScattering = new G4MultipleScattering();
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
      G4VProcess* aMultipleScattering = new G4MultipleScattering();
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
     G4VProcess* aMultipleScattering = new G4MultipleScattering();
     G4VProcess* anIonisation        = new G4hIonisation();
     //
     // add processes
     pmanager->AddProcess(anIonisation);
     pmanager->AddProcess(aMultipleScattering);
     //
     // set ordering for AlongStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,1);
     pmanager->SetProcessOrdering(anIonisation,        idxAlongStep,2);
     //
     // set ordering for PostStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep,1);
     pmanager->SetProcessOrdering(anIonisation,        idxPostStep,2);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4Decay.hh"

void ExN03PhysicsList::ConstructGeneral()
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

void ExN03PhysicsList::SetCuts()
{
  // reactualise cutValues
  if (currentDefaultCut != defaultCutValue)
    {
     if(cutForGamma    == currentDefaultCut) cutForGamma    = defaultCutValue;
     if(cutForElectron == currentDefaultCut) cutForElectron = defaultCutValue;
     if(cutForProton   == currentDefaultCut) cutForProton   = defaultCutValue;
     currentDefaultCut = defaultCutValue;
    }
    
  if (verboseLevel >0){
    G4cout << "ExN03PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;    
  }  

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForElectron, "e+");
 
  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton 
  SetCutValue(cutForProton, "proton");
  SetCutValue(cutForProton, "anti_proton");
  
  SetCutValueForOthers(defaultCutValue);

  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void ExN03PhysicsList::SetCutForGamma(G4double cut)
{
  ResetCuts();
  cutForGamma = cut;
}

void ExN03PhysicsList::SetCutForElectron(G4double cut)
{
  ResetCuts();
  cutForElectron = cut;
}

void ExN03PhysicsList::SetCutForProton(G4double cut)
{
  ResetCuts();
  cutForProton = cut;
}

G4double ExN03PhysicsList::GetCutForGamma() const
{
  return cutForGamma;
}

G4double ExN03PhysicsList::GetCutForElectron() const
{
  return cutForElectron;
}

G4double ExN03PhysicsList::GetCutForProton() const
{
  return cutForGamma;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....











