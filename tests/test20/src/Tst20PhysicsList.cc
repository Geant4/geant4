// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst20PhysicsList.cc,v 1.3 2001-05-25 13:19:47 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------

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

//#include "G4EnergyLossTables.hh"
//#include "G4Material.hh"
//#include "G4RunManager.hh"

//#include "G4UImanager.hh"

Tst20PhysicsList::Tst20PhysicsList(): G4VUserPhysicsList()
{
  defaultCutValue = 0.1*mm;
  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
  SetVerboseLevel(1);
  physicsListMessenger = new Tst20PhysicsListMessenger(this);
}


Tst20PhysicsList::~Tst20PhysicsList()
{
  delete physicsListMessenger;
}


void Tst20PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
}


void Tst20PhysicsList::ConstructBosons()
{
  
  // gamma
  G4Gamma::GammaDefinition();
  
}

void Tst20PhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

void Tst20PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
}


void Tst20PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* processManager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {

      // gamma, polarised processes
      processManager->AddDiscreteProcess(new G4LowEnergyPolarizedCompton);
      processManager->AddDiscreteProcess(new G4LowEnergyGammaConversion);

      lowEPhotoelProcess = new G4LowEnergyPhotoElectric();
      processManager->AddDiscreteProcess(lowEPhotoelProcess);

      processManager->AddDiscreteProcess(new G4LowEnergyRayleigh);
      
    } else if (particleName == "e-") {
      //electron
      processManager->AddProcess(new G4MultipleScattering,-1, 1,1);

      lowEIoniProcess = new G4LowEnergyIonisation();
      processManager->AddProcess(lowEIoniProcess, -1,  2, 2);

      lowEBremProcess = new G4LowEnergyBremsstrahlung();
      processManager->AddProcess(lowEBremProcess, -1, -1, 3);

    } else if (particleName == "e+") {
      //positron
      processManager->AddProcess(new G4MultipleScattering,-1, 1,1);
      processManager->AddProcess(new G4eIonisation,      -1, 2,2);
      processManager->AddProcess(new G4eBremsstrahlung,   -1,-1,3);
      processManager->AddProcess(new G4eplusAnnihilation,  0,-1,4);
      
    } 
  }
}


void Tst20PhysicsList::SetGELowLimit(G4double lowCut)
{
  if (verboseLevel >0){
    G4cout << "Tst20PhysicsList::SetCuts:";
    G4cout << "Gamma and Electron cut in energy: " 
	   << lowCut*MeV << " (MeV)" << G4endl;
  }  

  G4Gamma::SetEnergyRange(lowCut,1e5);
  G4Electron::SetEnergyRange(lowCut,1e5);
  G4Positron::SetEnergyRange(lowCut,1e5);

}

void Tst20PhysicsList::SetGammaLowLimit(G4double lowCut)
{
  if (verboseLevel >0){
    G4cout << "Tst20PhysicsList::SetCuts:";
    G4cout << "Gamma cut in energy: " 
	   << lowCut*MeV << " (MeV)" << G4endl;
  }  

  G4Gamma::SetEnergyRange(lowCut,1e5);

}

void Tst20PhysicsList::SetElectronLowLimit(G4double lowCut)
{
  if (verboseLevel >0){

    G4cout << "Tst20PhysicsList::SetCuts:";
    G4cout << "Electron cut in energy: " << lowCut*MeV << " (MeV)" << G4endl;

  }  

  G4Electron::SetEnergyRange(lowCut,1e5);

}
void Tst20PhysicsList::SetGammaCut(G4double val)
{
  ResetCuts();
  cutForGamma = val;
}


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


void Tst20PhysicsList::SetLowEnSecPhotCut(G4double cut){
  
  G4cout << "Low energy secondary photons cut is now set to: "
	 << cut*MeV << " (MeV)"<< G4endl;
  G4cout << "for LowEnergyBremsstrahlung, LowEnergyPhotoElectric, LowEnergyIonisation"
	 <<G4endl;
  lowEBremProcess->SetCutForLowEnSecPhotons(cut);
  lowEPhotoelProcess->SetCutForLowEnSecPhotons(cut);
  lowEIoniProcess->SetCutForLowEnSecPhotons(cut);
}

void Tst20PhysicsList::SetLowEnSecElecCut(G4double cut){
  
  G4cout << "Low energy secondary electrons cut is now set to: "
	 << cut*MeV << " (MeV)" << G4endl;
  //  G4cout <<"for processes LowEnergyBremsstrahlung, LowEnergyPhotoElectric, LowEnergyIonisation"<<G4endl;
  G4cout <<"for LowEnergyIonisation"<<G4endl;
  //  G4LowEnergyPhotoElectric::SetCutForLowEnSecElectrons(cut);
  lowEIoniProcess->SetCutForLowEnSecElectrons(cut);
}
