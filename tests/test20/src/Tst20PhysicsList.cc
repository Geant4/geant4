// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20PhysicsList.cc,v 1.2 2001-05-25 12:50:07 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Tst20PhysicsList.hh"
#include "Tst20PhysicsListMessenger.hh"
#include "Tst20DetectorConstruction.hh"
#include "G4LowEnergyPolarizedCompton.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"

//#include "G4EnergyLossTables.hh"
//#include "G4Material.hh"
//#include "G4RunManager.hh"

//#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20PhysicsList::Tst20PhysicsList()
: G4VUserPhysicsList()
{
  defaultCutValue = 0.1*mm;
  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
  SetVerboseLevel(1);
  physicsListMessenger = new Tst20PhysicsListMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Tst20PhysicsList::~Tst20PhysicsList()
{
  delete physicsListMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20PhysicsList::ConstructBosons()
{
  
  // gamma
  G4Gamma::GammaDefinition();
  
}
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20PhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4LowEnergyCompton.hh"
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyRayleigh.hh"

// e+
#include "G4MultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {

      // gamma         
      //      pmanager->AddDiscreteProcess(new G4LowEnergyCompton);
      pmanager->AddDiscreteProcess(new G4LowEnergyPolarizedCompton);
      pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion);

      LePeprocess = new G4LowEnergyPhotoElectric();
      pmanager->AddDiscreteProcess(LePeprocess);

      pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh);
      
    } else if (particleName == "e-") {
      //electron
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);

      LeIoprocess = new G4LowEnergyIonisation();
      pmanager->AddProcess(LeIoprocess, -1,  2, 2);

      LeBrprocess = new G4LowEnergyBremsstrahlung();
      pmanager->AddProcess(LeBrprocess, -1, -1, 3);

    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4eIonisation,      -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,   -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,  0,-1,4);
      
    } 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20PhysicsList::SetGELowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "Tst20PhysicsList::SetCuts:";
    G4cout << "Gamma and Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  G4Gamma::SetEnergyRange(lowcut,1e5);
  G4Electron::SetEnergyRange(lowcut,1e5);
  G4Positron::SetEnergyRange(lowcut,1e5);

}

void Tst20PhysicsList::SetGammaLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "Tst20PhysicsList::SetCuts:";
    G4cout << "Gamma cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  G4Gamma::SetEnergyRange(lowcut,1e5);

}

void Tst20PhysicsList::SetElectronLowLimit(G4double lowcut)
{
  if (verboseLevel >0){

    G4cout << "Tst20PhysicsList::SetCuts:";
    G4cout << "Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;

  }  

  G4Electron::SetEnergyRange(lowcut,1e5);

}
void Tst20PhysicsList::SetGammaCut(G4double val)
{
  ResetCuts();
  cutForGamma = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20PhysicsList::SetElectronCut(G4double val)
{
  ResetCuts();
  cutForElectron = val;
}

void Tst20PhysicsList::SetCuts(){

   SetCutValue(cutForGamma,"gamma");
   SetCutValue(cutForElectron,"e-");
   SetCutValue(cutForElectron,"e+");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Tst20PhysicsList::SetLowEnSecPhotCut(G4double cut){
  
  G4cout<<"Low energy secondary photons cut is now set to: "<<cut*MeV<<" (MeV)"<<G4endl;
  G4cout<<"for processes LowEnergyBremsstrahlung, LowEnergyPhotoElectric, LowEnergyIonisation"<<G4endl;
  LeBrprocess->SetCutForLowEnSecPhotons(cut);
  LePeprocess->SetCutForLowEnSecPhotons(cut);
  LeIoprocess->SetCutForLowEnSecPhotons(cut);
}

void Tst20PhysicsList::SetLowEnSecElecCut(G4double cut){
  
  G4cout<<"Low energy secondary electrons cut is now set to: "<<cut*MeV<<" (MeV)"<<G4endl;
  //  G4cout<<"for processes LowEnergyBremsstrahlung, LowEnergyPhotoElectric, LowEnergyIonisation"<<G4endl;
  G4cout<<"for processes LowEnergyIonisation"<<G4endl;
  //  G4LowEnergyPhotoElectric::SetCutForLowEnSecElectrons(cut);
  LeIoprocess->SetCutForLowEnSecElectrons(cut);
}













