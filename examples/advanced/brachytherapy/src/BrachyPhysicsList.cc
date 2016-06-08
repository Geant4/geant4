//    **********************************
//    *                                *
//    *     BrachyPhysicsList.cc       *
//    *                                *
//    **********************************

#include "BrachyPhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"              

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

BrachyPhysicsList::BrachyPhysicsList():  G4VUserPhysicsList()
{
  defaultCutValue = 1*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;
  
  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

BrachyPhysicsList::~BrachyPhysicsList()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BrachyPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BrachyPhysicsList::ConstructBosons()
{ 
  // gamma
  G4Gamma::GammaDefinition();

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BrachyPhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BrachyPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4MultipleScattering.hh"
// gamma
#include "G4LowEnergyRayleigh.hh" 
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"  
#include "G4LowEnergyGammaConversion.hh" 
// e-
#include "G4LowEnergyIonisation.hh" 
#include "G4LowEnergyBremsstrahlung.hh" 
// e+
#include "G4eIonisation.hh" 
#include "G4eBremsstrahlung.hh" 
#include "G4eplusAnnihilation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BrachyPhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    
    //processes
    lowePhot = new  G4LowEnergyPhotoElectric("LowEnPhotoElec");
    loweIon  = new G4LowEnergyIonisation("LowEnergyIoni");
    loweBrem = new G4LowEnergyBremsstrahlung("LowEnBrem");
    
    if (particleName == "gamma") {
      //gamma      
      pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh);
      pmanager->AddDiscreteProcess(lowePhot);
      pmanager->AddDiscreteProcess(new G4LowEnergyCompton);
      pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion);
      
    } else if (particleName == "e-") {
      //electron
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager->AddProcess(loweIon,     -1, 2,2);
      pmanager->AddProcess(loweBrem,    -1,-1,3);      
      
    } else if (particleName == "e+") {
      //positron      
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);      
      
    }
  }  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BrachyPhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "BrachyPhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }  
  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  
  SetCutValueForOthers(defaultCutValue);
  
  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BrachyPhysicsList::SetGammaLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "BrachyPhysicsList::SetCuts:";
    G4cout << "Gamma cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  
  
  G4Gamma::SetEnergyRange(lowcut,1e5);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BrachyPhysicsList::SetElectronLowLimit(G4double lowcut)
{
  if (verboseLevel >0){

    G4cout << "BrachyPhysicsList::SetCuts:";
    G4cout << "Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;

  }  

  G4Electron::SetEnergyRange(lowcut,1e5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BrachyPhysicsList::SetGELowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "BrachyPhysicsList::SetCuts:";
    G4cout << "Gamma and Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  G4Gamma::SetEnergyRange(lowcut,1e5);
  G4Electron::SetEnergyRange(lowcut,1e5);
  G4Positron::SetEnergyRange(lowcut,1e5);

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BrachyPhysicsList::SetGammaCut(G4double val)
{
  ResetCuts();
  cutForGamma = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BrachyPhysicsList::SetElectronCut(G4double val)
{
  //  ResetCuts();
  cutForElectron = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BrachyPhysicsList::SetPositronCut(G4double val)
{
  //  ResetCuts();
  cutForPositron = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BrachyPhysicsList::SetLowEnSecPhotCut(G4double cut){
  
  G4cout<<"Low energy secondary photons cut is now set to: "<<cut*MeV<<" (MeV)"<<G4endl;
  G4cout<<"for processes LowEnergyPhotoElectric, LowEnergyBremsstrahlung, LowEnergyIonisation"<<G4endl;
  lowePhot->SetCutForLowEnSecPhotons(cut);
  loweIon->SetCutForLowEnSecPhotons(cut);
  loweBrem->SetCutForLowEnSecPhotons(cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void BrachyPhysicsList::SetLowEnSecElecCut(G4double cut){
  
  G4cout<<"Low energy secondary electrons cut is now set to: "<<cut*MeV<<" (MeV)"<<G4endl;
  G4cout<<"for processes LowEnergyIonisation"<<G4endl;
  loweIon->SetCutForLowEnSecElectrons(cut);

}

