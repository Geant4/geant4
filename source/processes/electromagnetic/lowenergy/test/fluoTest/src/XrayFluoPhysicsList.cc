//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: XrayFluoPhysicsList.cc
// GEANT4 tag $Name: 
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "XrayFluoPhysicsList.hh"
#include "XrayFluoPhysicsListMessenger.hh"
#include "XrayFluoDetectorConstruction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ios.hh"
#include "XrayFluoDetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPhysicsList::XrayFluoPhysicsList(XrayFluoDetectorConstruction* p)
: G4VUserPhysicsList()
{
  pDet = p;
  defaultCutValue = 0.00001*mm;

  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
  cutForProton    = 0.00001*mm;
  SetVerboseLevel(1);
  physicsListMessenger = new XrayFluoPhysicsListMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPhysicsList::~XrayFluoPhysicsList()
{
  delete physicsListMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructBarions();
  ConstructIons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPhysicsList::ConstructBosons()
{
  
  // gamma
  G4Gamma::GammaDefinition();
  
}
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPhysicsList::ConstructBarions()
{
  G4Proton::ProtonDefinition();
}
void XrayFluoPhysicsList::ConstructIons()
{
//  Ions
 G4Alpha::AlphaDefinition();
}

void XrayFluoPhysicsList::ConstructProcess()
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
#include "G4hLowEnergyIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {

      // gamma         
      pmanager->AddDiscreteProcess(new G4LowEnergyCompton);
     
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
    else if (particleName == "proton") {
      //proton
      pmanager->AddProcess(new G4MultipleScattering,-1,1,1);
      pmanager->AddProcess(new G4hLowEnergyIonisation,-1, 2,2);
    }
    else if (   particleName == "alpha" )
      {
	
	pmanager->AddProcess(new G4MultipleScattering,-1,1,1);
	G4hLowEnergyIonisation* iIon = new G4hLowEnergyIonisation() ;
	pmanager->AddProcess(iIon,-1,2,2);
      }
   
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPhysicsList::SetGELowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "XrayFluoPhysicsList::SetCuts:";
    G4cout << "Gamma and Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  G4Gamma::SetEnergyRange(lowcut,1e5);
  G4Electron::SetEnergyRange(lowcut,1e5);
  G4Positron::SetEnergyRange(lowcut,1e5);

}

void XrayFluoPhysicsList::SetGammaLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "XrayFluoPhysicsList::SetCuts:";
    G4cout << "Gamma cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  G4Gamma::SetEnergyRange(lowcut,1e5);

}

void XrayFluoPhysicsList::SetElectronLowLimit(G4double lowcut)
{
  if (verboseLevel >0){

    G4cout << "XrayFluoPhysicsList::SetCuts:";
    G4cout << "Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;

  }  

  G4Electron::SetEnergyRange(lowcut,1e5);

}
void XrayFluoPhysicsList::SetGammaCut(G4double val)
{
  ResetCuts();
  cutForGamma = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPhysicsList::SetElectronCut(G4double val)
{
  ResetCuts();
  cutForElectron = val;
}

void XrayFluoPhysicsList::SetCuts(){

   SetCutValue(cutForGamma,"gamma");
   SetCutValue(cutForElectron,"e-");
   SetCutValue(cutForElectron,"e+");
   SetCutValue(cutForProton, "proton");
   SetCutValueForOthers(cutForProton);
   if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPhysicsList::SetLowEnSecPhotCut(G4double cut){
  
  G4cout<<"Low energy secondary photons cut is now set to: "<<cut*MeV<<" (MeV)"<<G4endl;
  G4cout<<"for processes LowEnergyBremsstrahlung, LowEnergyPhotoElectric, LowEnergyIonisation"<<G4endl;
  LeBrprocess->SetCutForLowEnSecPhotons(cut);
  LePeprocess->SetCutForLowEnSecPhotons(cut);
  LeIoprocess->SetCutForLowEnSecPhotons(cut);
}

void XrayFluoPhysicsList::SetLowEnSecElecCut(G4double cut){
  
  G4cout<<"Low energy secondary electrons cut is now set to: "<<cut*MeV<<" (MeV)"<<G4endl;
 
  G4cout<<"for processes LowEnergyIonisation"<<G4endl;
 
  LeIoprocess->SetCutForLowEnSecElectrons(cut);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPhysicsList::SetProtonCut(G4double val)
{
  //  ResetCuts();
  cutForProton = val;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XrayFluoPhysicsList::SetCutsByEnergy(G4double val)
{
  G4ParticleTable* theXrayFluoParticleTable =  G4ParticleTable::GetParticleTable();
  G4Material* currMat = pDet->GetSampleMaterial();
  G4ParticleDefinition* part;
  G4double cut;

  part = theXrayFluoParticleTable->FindParticle("e-");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "e-");
 
  part = theXrayFluoParticleTable->FindParticle("proton");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "proton");
 
}




































