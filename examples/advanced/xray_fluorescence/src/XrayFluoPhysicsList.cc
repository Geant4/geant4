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
// $Id: XrayFluoPhysicsList.cc
// GEANT4 tag $Name: xray_fluo-V03-02-00
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 28 Nov 2001 Elena Guardincerri     Created
//
// -------------------------------------------------------------------

#include "G4ParticleDefinition.hh" 
#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"
#include "XrayFluoPhysicsList.hh"
#include "XrayFluoPhysicsListMessenger.hh"
#include "XrayFluoDetectorConstruction.hh"
#include "XrayFluoPlaneDetectorConstruction.hh"
#include "XrayFluoMercuryDetectorConstruction.hh"

/////////////////////////////////////////
#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"
//#include "G4MesonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4BaryonConstructor.hh"
/////////////////////////////////////////

#include "G4StepLimiter.hh"

//#include "G4Region.hh"
//#include "G4RegionStore.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPhysicsList::XrayFluoPhysicsList(XrayFluoDetectorConstruction* p)
  : G4VUserPhysicsList(),pDet(0),planeDet(0),mercuryDet(0)

{
  pDet = p;

  SetGELowLimit(250*eV);

  defaultCutValue = 10e-6*mm;

  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
  cutForProton    = 0.00001*mm;
  SetVerboseLevel(1);
  physicsListMessenger = new XrayFluoPhysicsListMessenger(this);

  /*
  G4String regName = "SampleRegion";
  G4double cutValue = 0.000001 * mm;  
  G4Region* reg = G4RegionStore::GetInstance()->GetRegion(regName);

  G4ProductionCuts* cuts = new G4ProductionCuts;
  cuts->SetProductionCut(cutValue);
  reg->SetProductionCuts(cuts);
  */


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPhysicsList::XrayFluoPhysicsList(XrayFluoPlaneDetectorConstruction* p)
  : G4VUserPhysicsList(),pDet(0),planeDet(0),mercuryDet(0)
{
  planeDet = p;

  //  SetGELowLimit(250*eV);

  defaultCutValue = 0.000001*mm;

  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
  cutForProton    = 0.001*mm;
  SetVerboseLevel(1);
  physicsListMessenger = new XrayFluoPhysicsListMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

XrayFluoPhysicsList::XrayFluoPhysicsList(XrayFluoMercuryDetectorConstruction* p)
: G4VUserPhysicsList(),pDet(0),planeDet(0),mercuryDet(0)
{
  mercuryDet = p;

  //  SetGELowLimit(250*eV);

  defaultCutValue = 0.000001*mm;

  cutForGamma = defaultCutValue;
  cutForElectron = defaultCutValue;
  cutForProton    = 0.001*mm;
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
  G4BaryonConstructor baryon;
  baryon.ConstructParticle();
  G4Proton::ProtonDefinition();
}
void XrayFluoPhysicsList::ConstructIons()
{
//  Ions
  G4IonConstructor ions;
  ions.ConstructParticle();
  //  G4Alpha::AlphaDefinition();
}

void XrayFluoPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


// gamma from standard

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4RayleighScattering.hh" 

// gamma from Lowenergy

#include "G4LivermorePhotoElectricModel.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4LivermorePolarizedComptonModel.hh" // alternative for polarized photons
#include "G4LivermoreGammaConversionModel.hh"
#include "G4LivermoreRayleighModel.hh"
#include "G4LivermorePolarizedRayleighModel.hh" // alternative for polarized photons


// e+ - e- from standard

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4UniversalFluctuation.hh"
#include "G4eplusAnnihilation.hh"

// e+ - e- from Lowenergy

#include "G4LivermoreIonisationModel.hh"
#include "G4LivermoreBremsstrahlungModel.hh"


// protons from standard

#include "G4hLowEnergyIonisation.hh"
#include "G4hMultipleScattering.hh"

// options

#include "G4LossTableManager.hh"
#include "G4EmProcessOptions.hh"

//#include "G4hIonisation.hh"


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
      
      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      G4LivermoreComptonModel* theLivermoreComptonModel = new G4LivermoreComptonModel();
      theComptonScattering->SetModel(theLivermoreComptonModel);
      pmanager->AddDiscreteProcess(theComptonScattering);        
      
      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      theLivermorePhotoElectricModel = new G4LivermorePhotoElectricModel();
      theLivermorePhotoElectricModel->ActivateAuger(true);
      thePhotoElectricEffect->SetModel(theLivermorePhotoElectricModel);
      pmanager->AddDiscreteProcess(thePhotoElectricEffect);
     

      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      G4LivermoreGammaConversionModel* theLivermoreGammaConversionModel = new G4LivermoreGammaConversionModel();
      theGammaConversion->SetModel(theLivermoreGammaConversionModel);
      pmanager->AddDiscreteProcess(theGammaConversion);

      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      G4LivermoreRayleighModel* theRayleighModel = new G4LivermoreRayleighModel();
      theRayleigh->SetModel(theRayleighModel);
      pmanager->AddDiscreteProcess(theRayleigh);
      
      pmanager->AddProcess(new G4StepLimiter(), -1, -1, 5);


    } else if (particleName == "e-") {
      //electron 
            
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc, -1, 1, 1);
      
      G4eIonisation* eIoni = new G4eIonisation();
      G4LivermoreIonisationModel* theLivermoreIonisationModel = new G4LivermoreIonisationModel();
      theLivermoreIonisationModel->ActivateAuger(true);
      //      theLivermoreIonisationModel->SetUseAtomicDeexcitation(true);
      eIoni->SetEmModel(theLivermoreIonisationModel);
//      eIoni->SetStepFunction(0.2, 100*um); //     
      pmanager->AddProcess(eIoni, -1, 2, 2);
      
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
      eBrem->SetEmModel(new G4LivermoreBremsstrahlungModel());
      pmanager->AddProcess(eBrem, -1,-3, 3);
      pmanager->AddProcess(new G4StepLimiter(), -1, -1, 4);
      
      // Ionisation


      } else if (particleName == "e+") {
      //positron
      
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc, -1, 1, 1);
      
      pmanager->AddProcess(new G4eIonisation(), -1, 2, 2);



      pmanager->AddProcess(new G4eBremsstrahlung,   -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,  0,-1,4);
      
    }  
    else if (particleName == "proton") {
      //proton
      G4hLowEnergyIonisation* hIoni = new G4hLowEnergyIonisation();
      hIoni->SetFluorescence(true);
      hIoni->SelectShellIonisationCS("analytical");
      pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);
      pmanager->AddProcess(hIoni,-1, 2,2);
    }
    else if (   particleName == "alpha" )
      {
	
	pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);
	G4hLowEnergyIonisation* hIoni = new G4hLowEnergyIonisation();
	hIoni->SetFluorescence(true);	
	pmanager->AddProcess(hIoni,-1,2,2);
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

  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowcut,1e5);


}

void XrayFluoPhysicsList::SetGammaLowLimit(G4double lowcut)
{
  if (verboseLevel >0){
    G4cout << "XrayFluoPhysicsList::SetCuts:";
    G4cout << "Gamma cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;
  }  

  SetGELowLimit(lowcut);

//   G4Gamma::SetEnergyRange(lowcut,1e5);

}

void XrayFluoPhysicsList::SetElectronLowLimit(G4double lowcut)
{
  if (verboseLevel >0){

    G4cout << "XrayFluoPhysicsList::SetCuts:";
    G4cout << "Electron cut in energy: " << lowcut*MeV << " (MeV)" << G4endl;

  }  

  SetGELowLimit(lowcut);

//   G4Electron::SetEnergyRange(lowcut,1e5);

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

  G4double lowlimit=20*eV;
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit, 100.*GeV);
  SetCutValue(cutForGamma,"gamma");
  SetCutValue(cutForElectron,"e-");
  SetCutValue(cutForElectron,"e+");
  SetCutValue(cutForProton, "proton");
  //SetCutValueForOthers(cutForProton);
  if (verboseLevel>0) DumpCutValuesTable();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPhysicsList::SetLowEnSecPhotCut(G4double cut){

  cut=0;
  
  G4cout<<"Low energy secondary photons cut is now set to: "<<cut/MeV<<" (MeV)"<<G4endl;
  G4cout<<"for processes LowEnergyBremsstrahlung, LowEnergyPhotoElectric, LowEnergyIonisation"<<G4endl;
//  LeBrprocess->SetCutForLowEnSecPhotons(cut);
//  theLivermorePhotoElectricModel->SetCutForLowEnSecPhotons(cut);
//  LeIoprocess->SetCutForLowEnSecPhotons(cut);

}

/*
void XrayFluoPhysicsList::SetLowEnSecElecCut(G4double cut){
  
  G4cout<<"Low energy secondary electrons cut is now set to: "<<cut/MeV<<" (MeV)"<<G4endl;
  
  G4cout<<"for processes LowEnergyIonisation"<<G4endl;
  
  LeIoprocess->SetCutForLowEnSecElectrons(cut);

}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void XrayFluoPhysicsList::SetProtonCut(G4double val)
{
  ResetCuts();
  cutForProton = val;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


void XrayFluoPhysicsList::SetCutsByEnergy(G4double val)
{
  G4ParticleTable* theXrayFluoParticleTable =  G4ParticleTable::GetParticleTable();

  G4Material* currMat = 0;
    
  if(pDet){
    currMat = pDet->XrayFluoDetectorConstruction::GetSampleMaterial();
  }
  else if(planeDet){
    currMat = planeDet->XrayFluoPlaneDetectorConstruction::GetPlaneMaterial(); 
  }
  
  else if(mercuryDet){
    currMat = mercuryDet->GetMercuryMaterial();
  }

  G4ParticleDefinition* part;
  G4double cut;
  
  part = theXrayFluoParticleTable->FindParticle("e-");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "e-");
  
  part = theXrayFluoParticleTable->FindParticle("proton");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "proton");
 
  part = theXrayFluoParticleTable->FindParticle("gamma");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "gamma");

}
























