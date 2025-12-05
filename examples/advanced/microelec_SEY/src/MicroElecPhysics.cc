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
// -------------------------------------------------------------------
// -------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "MicroElecPhysics.hh"
#include "G4SystemOfUnits.hh"


// Geant4-MicroElec MODELS

#include "G4MicroElecElastic.hh"
#include "G4MicroElecElasticModel_new.hh"

#include "G4MicroElecInelastic.hh"
#include "G4MicroElecInelasticModel_new.hh"

#include "G4MicroElecLOPhononScattering.hh"
#include "G4MicroElecLOPhononModel.hh"
#include "G4MicroElecSurface.hh"

//

#include "G4LossTableManager.hh"
#include "G4EmConfigurator.hh"
#include "G4VEmModel.hh"
#include "G4DummyModel.hh"
#include "G4eIonisation.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4BraggModel.hh"
#include "G4BraggIonModel.hh"
#include "G4BetheBlochModel.hh"
#include "G4UrbanMscModel.hh"
#include "G4MollerBhabhaModel.hh"
#include "G4IonFluctuations.hh"
#include "G4UniversalFluctuation.hh"

#include "G4MicroElecCapture.hh"

#include "G4UAtomicDeexcitation.hh"



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicroElecPhysics::MicroElecPhysics():  G4VUserPhysicsList()
{
  defaultCutValue = 1*micrometer;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;
  cutForPositron  = defaultCutValue;
  cutForProton    = defaultCutValue;
  
  SetVerboseLevel(1);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

MicroElecPhysics::~MicroElecPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicroElecPhysics::ConstructParticle()
{
  ConstructBosons();
  ConstructLeptons();
  ConstructBarions();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicroElecPhysics::ConstructBosons()
{ 
  // gamma
  G4Gamma::GammaDefinition();
}
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicroElecPhysics::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicroElecPhysics::ConstructBarions()
{
  //  baryons
  G4Proton::ProtonDefinition();
  G4GenericIon::GenericIonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicroElecPhysics::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicroElecPhysics::ConstructEM()
{

  G4EmParameters* param = G4EmParameters::Instance();
  param->SetBuildCSDARange(true);
  param->SetMscStepLimitType(fUseSafety);
  param->RegionsMicroElec();
  // physicList ISS
  param->SetDefaults();
  param->SetMinEnergy(0.1*eV);
  param->SetMaxEnergy(10 * TeV);
  param->SetLowestElectronEnergy(0 * eV); //<--- Energy cut in the vacuum!!! =0eV to get right SEY
  param->SetNumberOfBinsPerDecade(20);
  param->ActivateAngularGeneratorForIonisation(true);
  param->SetAugerCascade(true);//*/
  
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  
  while( (*particleIterator)() )
  {

    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    // *********************************
    // 1) Processes for the World region
    // *********************************

    if (particleName == "e-") {

      // STANDARD msc is active in the world
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->AddEmModel(1, new G4UrbanMscModel());
      pmanager->AddProcess(msc, -1, 1, -1);

      // STANDARD ionisation is active in the world
      G4eIonisation* eion = new G4eIonisation();
      pmanager->AddProcess(eion, -1, 2, 2);

      // MicroElec elastic is not active in the world 
	  G4MicroElecElastic* theMicroElecElasticProcess = new G4MicroElecElastic("e-_G4MicroElecElastic");
      theMicroElecElasticProcess->SetEmModel(new G4DummyModel(),1);
      pmanager->AddDiscreteProcess(theMicroElecElasticProcess);
	    	    
      // MicroElec ionisation is not active in the world 
      G4MicroElecInelastic* microelecioni = new G4MicroElecInelastic("e-_G4Dielectrics");
      microelecioni->SetEmModel(new G4DummyModel(),1); 
      pmanager->AddDiscreteProcess(microelecioni);
	  
	  //Phonons for SiO2, AL2O3, BN
	  G4MicroElecLOPhononScattering* opticalPhonon = new G4MicroElecLOPhononScattering("e-_G4LOPhononScattering");
	  opticalPhonon->SetEmModel(new G4DummyModel(), 1);
	  pmanager->AddDiscreteProcess(opticalPhonon);


	  G4MicroElecSurface* MicroElecSurf = new G4MicroElecSurface("e-_G4MicroElecSurface");
	  MicroElecSurf->SetProcessManager(pmanager);
	  pmanager->AddDiscreteProcess(MicroElecSurf);//*/

      // Modified for capture in Al2O3 and SiO2 (Polaronic trapping effect)
      G4MicroElecCapture* ecap = new G4MicroElecCapture("Target",0.9*eV); 
	  pmanager->AddDiscreteProcess(ecap);//*/
     	    
    } else if ( particleName == "proton" ) {
		
      // multiple scatering for protons
      G4hMultipleScattering* msc = new G4hMultipleScattering();
      pmanager->AddProcess(msc, -1, 1, 1); 

      // STANDARD ionisation is active in the world 
      G4hIonisation* hion = new G4hIonisation();
      pmanager->AddProcess(hion, -1, 2, 2);

    } else if(particleName == "alpha") {
		
	  // multiple scatering for protons
      G4hMultipleScattering* msc = new G4hMultipleScattering();
      pmanager->AddProcess(msc, -1, 1, 1); 

      // STANDARD ionisation is active in the world
      G4ionIonisation* hion = new G4ionIonisation();
      pmanager->AddProcess(hion, -1, 2, 2);
      
    } else if (particleName == "GenericIon") {

	  // multiple scatering for protons
      G4hMultipleScattering* msc = new G4hMultipleScattering();
      pmanager->AddProcess(msc, -1, 1, 1); 
	  
      // STANDARD ionisation is active in the world 
      G4ionIonisation* hion = new G4ionIonisation();
      pmanager->AddProcess(hion, -1, 2, 2);

    } 
  }

  // **************************************
  // 2) Define processes for Target region 
  // **************************************

  // STANDARD EM processes should be inactivated when corresponding MicroElec processes are used
  // - STANDARD EM e- processes are inactivated below 100 MeV
  // - STANDARD EM proton & ion processes are inactivated below standEnergyLimit
  //
  G4EmConfigurator* em_config = G4LossTableManager::Instance()->EmConfigurator();

  G4VEmModel* mod;
  
  // *** e- ---------------------------------------------------------- 
  // ---> STANDARD EM processes are inactivated below 100 MeV
  G4UrbanMscModel* msc =  new G4UrbanMscModel();
  msc->SetActivationLowEnergyLimit(100*MeV);
  em_config->SetExtraEmModel("e-","msc",msc,"Target");
  
  mod = new G4MollerBhabhaModel();
  mod->SetActivationLowEnergyLimit(10*MeV);
  em_config->SetExtraEmModel("e-","eIoni",mod,"Target",0.0,10*TeV, new G4UniversalFluctuation());
  
  // ---> MicroElec processes activated
  mod = new G4MicroElecElasticModel_new();
  em_config->SetExtraEmModel("e-","e-_G4MicroElecElastic",mod,"Target",0.1*eV,100*MeV);
 
  mod = new G4MicroElecInelasticModel_new();
  em_config->SetExtraEmModel("e-","e-_G4Dielectrics",mod,"Target",0.1*eV,10*MeV);

  //Phonons LO pour sio2,  al2o3 and BN
  mod = new G4MicroElecLOPhononModel();
  em_config->SetExtraEmModel("e-", "e-_G4LOPhononScattering", mod, "Target", 0.1 * eV, 10 * MeV);//*/

  // *** proton ---------------------------------------------------------- 
  mod = new G4BetheBlochModel();
  mod->SetActivationLowEnergyLimit(10*MeV);
  em_config->SetExtraEmModel("proton","hIoni",mod,"Target",2*MeV,10*TeV, new G4IonFluctuations());
			     
  // ---> Dielectric processes activated
  mod = new G4MicroElecInelasticModel_new();
  mod->SetActivationLowEnergyLimit(100*eV);
  em_config->SetExtraEmModel("proton","p_G4Dielectrics",mod,"Target",100*eV,10*MeV);

  // *** alpha ----------------------------------------------------------
  mod = new G4BetheBlochModel();
  mod->SetActivationLowEnergyLimit(10*MeV);
  em_config->SetExtraEmModel("alpha","ionIoni",mod,"Target",10*MeV,10*TeV, new G4IonFluctuations());

  // *** ion ----------------------------------------------------------
  // ---> STANDARD EM processes inactivated below standEnergyLimit

  mod = new G4BetheBlochModel();
  mod->SetActivationLowEnergyLimit(10*MeV);
  em_config->SetExtraEmModel("GenericIon","ionIoni",mod,"Target",10*MeV,10*TeV, new G4IonFluctuations());
   
  // ---> Dielectric processes activated
  mod = new G4MicroElecInelasticModel_new();
  mod->SetActivationLowEnergyLimit(100*eV);
  em_config->SetExtraEmModel("GenericIon","ion_G4Dielectrics",mod,"Target",0.0,10*MeV);

  // Deexcitation
  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
  de->SetFluo(true);
  de->SetAuger(true);   
  de->SetPIXE(true);  
  de->InitialiseForNewRun();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void MicroElecPhysics::ConstructGeneral()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

