//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "FluoTestPhysicsList.hh"
#include "FluoTestPhysicsListMessenger.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4ios.hh"              

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestPhysicsList::FluoTestPhysicsList():  G4VUserPhysicsList()
{
  currentDefaultCut = defaultCutValue = 1.0*mm;
  cutForGamma       = defaultCutValue;
  cutForElectron    = defaultCutValue;
  cutForProton      = defaultCutValue;
 physicsListMessenger = new FluoTestPhysicsListMessenger(this);
 SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

FluoTestPhysicsList::~FluoTestPhysicsList()
{
 delete physicsListMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestPhysicsList::ConstructParticle()
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

void FluoTestPhysicsList::ConstructBosons()
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

void FluoTestPhysicsList::ConstructLeptons()
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

void FluoTestPhysicsList::ConstructMesons()
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

void FluoTestPhysicsList::ConstructBaryons()
{
//  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//#include "G4ComptonScattering.hh"
//#include "G4GammaConversion.hh"
//#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hIonisation.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyRayleigh.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void FluoTestPhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
  if (particleName == "gamma") {
    // gamma
     pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion());
    //pmanager->AddDiscreteProcess(new G4GammaConversion());
     pmanager->AddDiscreteProcess(new G4LowEnergyCompton());      
     //pmanager->AddDiscreteProcess(new G4ComptonScattering());    
     pmanager->AddDiscreteProcess(new G4LowEnergyPhotoElectric());
    // pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
  }else if (particleName == "e-") {

   
    //electron
      G4VProcess* theeminusMultipleScattering = new G4MultipleScattering();
      G4VProcess* theeminusIonisation         = new G4LowEnergyIonisation();
      //G4VProcess* theeminusIonisation         = new G4eIonisation();
      G4VProcess* theeminusBremsstrahlung     = new G4LowEnergyBremsstrahlung();
      //G4VProcess* theeminusBremsstrahlung     = new G4eBremsstrahlung();
       G4VProcess* theeminusRayleigh           = new G4LowEnergyRayleigh();
     
      //
      // add processes
      pmanager->AddProcess(theeminusMultipleScattering);
      pmanager->AddProcess(theeminusIonisation);
      pmanager->AddProcess(theeminusBremsstrahlung);
      pmanager->AddProcess(theeminusRayleigh);
      //      
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxAlongStep,1);
      pmanager->SetProcessOrdering(theeminusIonisation,         idxAlongStep,2);
      //
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxPostStep,1);
      pmanager->SetProcessOrdering(theeminusIonisation,         idxPostStep,2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung,     idxPostStep,3);
      pmanager->SetProcessOrdering(theeminusRayleigh,           idxPostStep,4);
    
    } else if (particleName == "e+") {
    //positron
      G4VProcess* theeplusMultipleScattering = new G4MultipleScattering();
      G4VProcess* theeplusIonisation         = new G4eIonisation();
      G4VProcess* theeplusAnnihilation       = new G4eplusAnnihilation();
      G4VProcess* theeplusBremsstrahlung     = new G4eBremsstrahlung();
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

void FluoTestPhysicsList::ConstructGeneral()
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

void FluoTestPhysicsList::SetGELowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "FluoTestPhysicsList::SetCuts:";
    G4cout << "Gamma and Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  G4Gamma::SetEnergyRange(lowcut,1e5);
  G4Electron::SetEnergyRange(lowcut,1e5);
  G4Positron::SetEnergyRange(lowcut,1e5);

}


void FluoTestPhysicsList::SetGammaLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "FluoTestPhysicsList::SetCuts:";
    G4cout << "Gamma cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  G4Gamma::SetEnergyRange(lowcut,1e5);

} 
void FluoTestPhysicsList::SetElectronLowLimit(G4double lowcut)
{
  if (verboseLevel >0){

    G4cout << "FluoTestPhysicsList::SetCuts:";
    G4cout << "Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;

  }  

  G4Electron::SetEnergyRange(lowcut,1e5);

} 

 void FluoTestPhysicsList::SetCuts()
{
 
  // reactualise cutValues
  if (currentDefaultCut != defaultCutValue)
    {
    
     if(cutForProton   == currentDefaultCut) cutForProton   = defaultCutValue;
     currentDefaultCut = defaultCutValue;
    }
    
  if (verboseLevel >0){
    G4cout << "FluoTestPhysicsList::SetCuts:";
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


void FluoTestPhysicsList::SetCutForGamma(G4double cut)
{
  ResetCuts();
  cutForGamma = cut;
}

void FluoTestPhysicsList::SetCutForElectron(G4double cut)
{
  ResetCuts();
  cutForElectron = cut;
}

void FluoTestPhysicsList::SetCutForProton(G4double cut)
{
  ResetCuts();
  cutForProton = cut;
}

G4double FluoTestPhysicsList::GetCutForGamma() const
{
  return cutForGamma;
}
G4double FluoTestPhysicsList::GetCutForElectron() const
{
  return cutForElectron;
}

G4double FluoTestPhysicsList::GetCutForProton() const
{
  return cutForProton;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....









