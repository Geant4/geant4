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
// add elastic scattering processes of proton, hydrogen, helium, alpha+, alpha++

#include "G4EmDNAPhysics.hh"

#include "G4SystemOfUnits.hh"

// *** Processes and models for Geant4-DNA
#include "G4DNAElectronSolvation.hh"
#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"
#include "G4DNAIonElasticModel.hh"

#include "G4DNAExcitation.hh"
#include "G4DNAAttachment.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DNAIonisation.hh"
#include "G4DNAChargeDecrease.hh"
#include "G4DNAChargeIncrease.hh"

// particles
#include "G4Electron.hh"
#include "G4Proton.hh"
#include "G4Alpha.hh"
#include "G4GenericIon.hh"
#include "G4DNAGenericIonsManager.hh"

// e+
#include "G4Positron.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

// gamma
#include "G4Gamma.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4KleinNishinaModel.hh"
#include "G4GammaConversion.hh"
#include "G4RayleighScattering.hh" 
#include "G4BetheHeitler5DModel.hh"

// utilities
#include "G4EmParameters.hh"
#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"
#include "G4EmBuilder.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysics::G4EmDNAPhysics(G4int ver, const G4String& name)
  : G4VPhysicsConstructor(name)
{
  // parameters for DNA and for option3 EM physics
  SetVerboseLevel(ver);
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetMinEnergy(10*CLHEP::eV);
  param->SetLowestElectronEnergy(0*CLHEP::eV);
  param->SetNumberOfBinsPerDecade(20);
  param->SetStepFunction(0.2, 100*CLHEP::um);
  param->SetMscStepLimitType(fUseDistanceToBoundary);
  param->SetLateralDisplacementAlg96(true);
  param->SetFluo(true);  
  param->SetAuger(true);  
  param->SetDeexcitationIgnoreCut(true);
  param->ActivateDNA();

  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysics::~G4EmDNAPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysics::ConstructParticle()
{
  // bosons
  G4Gamma::Gamma();

  // leptons
  G4Electron::Electron();
  G4Positron::Positron();
  
  // baryons
  G4Proton::Proton();

  G4GenericIon::GenericIonDefinition();
  G4Alpha::Alpha();

  G4DNAGenericIonsManager* genericIonsManager
    = G4DNAGenericIonsManager::Instance();
  //genericIonsManager->GetIon("alpha++");
  genericIonsManager->GetIon("alpha+");
  genericIonsManager->GetIon("helium");
  genericIonsManager->GetIon("hydrogen");
  //genericIonsManager->GetIon("carbon");
  //genericIonsManager->GetIon("nitrogen");
  //genericIonsManager->GetIon("oxygen");
  //genericIonsManager->GetIon("iron");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysics::ConstructProcess()
{
  if(verboseLevel > 1) {
    G4cout << "### " << GetPhysicsName() 
	   << " Construct Processes " << G4endl;
  }
  ConstructGammaPositronProcesses();

  G4PhysicsListHelper* helper = 
    G4PhysicsListHelper::GetPhysicsListHelper();
  G4DNAGenericIonsManager* genericIonsManager
    = G4DNAGenericIonsManager::Instance();

  // e-
  G4ParticleDefinition* part = G4Electron::Electron();

  // *** Solvation ***
  G4DNAElectronSolvation* pSolvation = 
    new G4DNAElectronSolvation("e-_G4DNAElectronSolvation");
  auto pSolvationModel = G4DNASolvationModelFactory::GetMacroDefinedModel();
  pSolvationModel->SetHighEnergyLimit(7.4*eV); // limit of the Champion's model
  pSolvation->SetEmModel(pSolvationModel);
  helper->RegisterProcess(pSolvation, part);
      
  // *** Elastic scattering ***
  G4DNAElastic* pElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
  pElasticProcess->SetEmModel(new G4DNAChampionElasticModel());   
  helper->RegisterProcess(pElasticProcess, part);

  // *** Excitation ***
  helper->RegisterProcess(new G4DNAExcitation("e-_G4DNAExcitation"), part);

  // *** Ionisation ***
  helper->RegisterProcess(new G4DNAIonisation("e-_G4DNAIonisation"), part);

  // *** Vibrational excitation ***
  helper->RegisterProcess(new G4DNAVibExcitation("e-_G4DNAVibExcitation"), part);
      
  // *** Attachment ***
  helper->RegisterProcess(new G4DNAAttachment("e-_G4DNAAttachment"), part); 

  // proton
  part = G4Proton::Proton();

  helper->RegisterProcess(new G4DNAElastic("proton_G4DNAElastic"), part);
  helper->RegisterProcess(new G4DNAExcitation("proton_G4DNAExcitation"), part);
  helper->RegisterProcess(new G4DNAIonisation("proton_G4DNAIonisation"), part);
  helper->RegisterProcess(new G4DNAChargeDecrease("proton_G4DNAChargeDecrease"), part);

  // hydrogen
  part = genericIonsManager->GetIon("hydrogen");

  helper->RegisterProcess(new G4DNAElastic("hydrogen_G4DNAElastic"), part);
  helper->RegisterProcess(new G4DNAExcitation("hydrogen_G4DNAExcitation"), part);
  helper->RegisterProcess(new G4DNAIonisation("hydrogen_G4DNAIonisation"), part);
  helper->RegisterProcess(new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease"), part);

  // alpha++
  part = G4Alpha::Alpha();

  helper->RegisterProcess(new G4DNAElastic("alpha_G4DNAElastic"), part);
  helper->RegisterProcess(new G4DNAExcitation("alpha_G4DNAExcitation"), part);
  helper->RegisterProcess(new G4DNAIonisation("alpha_G4DNAIonisation"), part);
  helper->RegisterProcess(new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease"), part);

  // alpha+
  part = genericIonsManager->GetIon("alpha+");

  helper->RegisterProcess(new G4DNAElastic("alpha+_G4DNAElastic"), part);
  helper->RegisterProcess(new G4DNAExcitation("alpha+_G4DNAExcitation"), part);
  helper->RegisterProcess(new G4DNAIonisation("alpha+_G4DNAIonisation"), part);
  helper->RegisterProcess(new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease"), part);
  helper->RegisterProcess(new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease"), part);

  // helium
  part = genericIonsManager->GetIon("helium");

  helper->RegisterProcess(new G4DNAElastic("helium_G4DNAElastic"), part);
  helper->RegisterProcess(new G4DNAExcitation("helium_G4DNAExcitation"), part);
  helper->RegisterProcess(new G4DNAIonisation("helium_G4DNAIonisation"), part);
  helper->RegisterProcess(new G4DNAChargeIncrease("helium_G4DNAChargeIncrease"), part);
    
  // other ions
  part = G4GenericIon::GenericIon();

  helper->RegisterProcess(new G4DNAIonisation("GenericIon_G4DNAIonisation"), part);
}

void G4EmDNAPhysics::ConstructGammaPositronProcesses()
{
  // this construction is based on G4EmStandardPhysics_option3
  G4EmBuilder::PrepareEMPhysics();

  G4PhysicsListHelper* helper = G4PhysicsListHelper::GetPhysicsListHelper();
  G4ParticleDefinition* part = G4Gamma::Gamma();

  // photoelectric effect - Livermore model 
  G4PhotoElectricEffect* thePEEffect = new G4PhotoElectricEffect();
  thePEEffect->SetEmModel(new G4LivermorePhotoElectricModel());
  helper->RegisterProcess(thePEEffect, part);

  // Compton scattering - Klein-Nishina
  G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
  theComptonScattering->SetEmModel(new G4KleinNishinaModel());
  helper->RegisterProcess(theComptonScattering, part);

  // gamma conversion - 5D model
  G4GammaConversion* theGammaConversion = new G4GammaConversion();
  theGammaConversion->SetEmModel(new G4BetheHeitler5DModel());
  helper->RegisterProcess(theGammaConversion, part);

  // Rayleigh scattering - Livermore model
  G4RayleighScattering* theRayleigh = new G4RayleighScattering();
  helper->RegisterProcess(theRayleigh, part);

  part = G4Positron::Positron();

  helper->RegisterProcess(new G4eMultipleScattering(), part);
  helper->RegisterProcess(new G4eIonisation(), part);
  helper->RegisterProcess(new G4eBremsstrahlung(), part);
  helper->RegisterProcess(new G4eplusAnnihilation(), part);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
