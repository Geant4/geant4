// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em1PhysicsList.cc,v 1.1 1999-10-11 13:07:46 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em1PhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"

#include "G4ParticleTable.hh"
#include "G4EnergyLossTables.hh"
#include "G4Material.hh"
#include "G4RunManager.hh"
#include "Em1DetectorConstruction.hh"
#include "Em1PhysicsListMessenger.hh"

#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1PhysicsList::Em1PhysicsList(Em1DetectorConstruction* p)
: G4VUserPhysicsList()
{
  pDet = p;

  defaultCutValue = 1.*mm;
  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
  cutForProton = defaultCutValue;

  SetVerboseLevel(1);
  pMes = new Em1PhysicsListMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em1PhysicsList::~Em1PhysicsList()
{ delete pMes; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBarions();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1PhysicsList::ConstructBosons()
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

void Em1PhysicsList::ConstructLeptons()
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

void Em1PhysicsList::ConstructMesons()
{
 //  mesons
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1PhysicsList::ConstructBarions()
{
//  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1PhysicsList::ConstructProcess()
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

void Em1PhysicsList::ConstructEM()
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
      pmanager->AddProcess(new G4eBremsstrahlung,   -1,-1,3);      

    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4eIonisation,       -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,   -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,  0,-1,4);
      
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon  
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4MuIonisation,      -1, 2,2);
      pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1,3);
      pmanager->AddProcess(new G4MuPairProduction,  -1,-1,4);       
     
    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddProcess(new G4MultipleScattering,-1,1,1);
      pmanager->AddProcess(new G4hIonisation,       -1,2,2);                  
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4Decay.hh"

void Em1PhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while ((*theParticleIterator)()){
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

void Em1PhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "Em1PhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue/mm << " (mm)" << endl;
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

void Em1PhysicsList::SetGammaCut(G4double val)
{
  ResetCuts();
  cutForGamma = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1PhysicsList::SetElectronCut(G4double val)
{
  ResetCuts();
  cutForElectron = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1PhysicsList::SetProtonCut(G4double val)
{
  ResetCuts();
  cutForProton = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1PhysicsList::SetCutsByEnergy(G4double val)
{
  G4ParticleTable* theParticleTable =  G4ParticleTable::GetParticleTable();
  G4Material* currMat = pDet->GetMaterial();
  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  G4ParticleDefinition* part;
  G4double cut;
  
  part = theParticleTable->FindParticle("e-");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "e-");

  part = theParticleTable->FindParticle("e+");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "e+");
  
  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton 
  part = theParticleTable->FindParticle("proton");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "proton");

  part = theParticleTable->FindParticle("anti_proton");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "anti_proton");
   
  SetCutValueForOthers(cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em1PhysicsList::GetRange(G4double val)
{
  G4ParticleTable* theParticleTable =  G4ParticleTable::GetParticleTable();
  G4Material* currMat = pDet->GetMaterial();

  G4ParticleDefinition* part;
  G4double cut;
  part = theParticleTable->FindParticle("e-");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  G4cout << "material : " << currMat->GetName() << endl;
  G4cout << "particle : " << part->GetParticleName() << endl;
  G4cout << "energy   : " << val / keV << " (keV)" << endl;
  G4cout << "range    : " << cut / mm << " (mm)" << endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


