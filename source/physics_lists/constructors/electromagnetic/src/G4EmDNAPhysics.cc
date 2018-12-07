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

#include "G4DNAGenericIonsManager.hh"

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
#include "G4GenericIon.hh"

// Warning : the following is needed in order to use EM Physics builders
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
#include "G4LivermoreComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"

#include "G4EmParameters.hh"
// end of warning

#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4PhysicsListHelper.hh"
#include "G4BuilderType.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAPhysics);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysics::G4EmDNAPhysics(G4int ver, const G4String&)
  : G4VPhysicsConstructor("G4EmDNAPhysics"), verbose(ver)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetFluo(true);  
  param->SetAuger(true);  
  param->SetAugerCascade(true);  
  param->SetDeexcitationIgnoreCut(true);

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

  G4DNAGenericIonsManager * genericIonsManager;
  genericIonsManager=G4DNAGenericIonsManager::Instance();
  genericIonsManager->GetIon("alpha++");
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
  if(verbose > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4PhysicsListHelper* pPhysicsHelper = G4PhysicsListHelper::GetPhysicsListHelper();

  auto pParticleIterator = GetParticleIterator();
  pParticleIterator->reset();
  while( (*pParticleIterator)() )
  {
    G4ParticleDefinition* pParticle = pParticleIterator->value();
    G4String particleName = pParticle->GetParticleName();

    if (particleName == "e-") {

      G4DNAElectronSolvation* pSolvation = new G4DNAElectronSolvation("e-_G4DNAElectronSolvation");
      auto pSolvationModel = G4DNASolvationModelFactory::GetMacroDefinedModel();
      pSolvationModel->SetHighEnergyLimit(7.4*eV); // limit of the Champion's model
      pSolvation->SetEmModel(pSolvationModel);
      pPhysicsHelper->RegisterProcess(pSolvation, pParticle);
      
      // *** Elastic scattering (two alternative models available) ***
      
      G4DNAElastic* pElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
      pElasticProcess->SetEmModel(new G4DNAChampionElasticModel());
      
      // or alternative model
      //theDNAElasticProcess->SetEmModel(new G4DNAScreenedRutherfordElasticModel());
      
      pPhysicsHelper->RegisterProcess(pElasticProcess, pParticle);

      // *** Excitation ***
      pPhysicsHelper->RegisterProcess(new G4DNAExcitation("e-_G4DNAExcitation"), pParticle);

      // *** Ionisation ***
      pPhysicsHelper->RegisterProcess(new G4DNAIonisation("e-_G4DNAIonisation"), pParticle);

      // *** Vibrational excitation ***
      pPhysicsHelper->RegisterProcess(new G4DNAVibExcitation("e-_G4DNAVibExcitation"), pParticle);
      
      // *** Attachment ***
      pPhysicsHelper->RegisterProcess(new G4DNAAttachment("e-_G4DNAAttachment"), pParticle); 
    
    } else if ( particleName == "proton" ) {
      pPhysicsHelper->RegisterProcess(new G4DNAElastic("proton_G4DNAElastic"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAExcitation("proton_G4DNAExcitation"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAIonisation("proton_G4DNAIonisation"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAChargeDecrease("proton_G4DNAChargeDecrease"), pParticle);

    } else if ( particleName == "hydrogen" ) {
      pPhysicsHelper->RegisterProcess(new G4DNAElastic("hydrogen_G4DNAElastic"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAExcitation("hydrogen_G4DNAExcitation"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAIonisation("hydrogen_G4DNAIonisation"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease"), pParticle);

    } else if ( particleName == "alpha" ) {
      pPhysicsHelper->RegisterProcess(new G4DNAElastic("alpha_G4DNAElastic"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAExcitation("alpha_G4DNAExcitation"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAIonisation("alpha_G4DNAIonisation"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease"), pParticle);

    } else if ( particleName == "alpha+" ) {
      pPhysicsHelper->RegisterProcess(new G4DNAElastic("alpha+_G4DNAElastic"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAExcitation("alpha+_G4DNAExcitation"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAIonisation("alpha+_G4DNAIonisation"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease"), pParticle);

    } else if ( particleName == "helium" ) {
      pPhysicsHelper->RegisterProcess(new G4DNAElastic("helium_G4DNAElastic"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAExcitation("helium_G4DNAExcitation"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAIonisation("helium_G4DNAIonisation"), pParticle);
      pPhysicsHelper->RegisterProcess(new G4DNAChargeIncrease("helium_G4DNAChargeIncrease"), pParticle);
    
    } else if ( particleName == "GenericIon" ) {
      pPhysicsHelper->RegisterProcess(new G4DNAIonisation("GenericIon_G4DNAIonisation"), pParticle);

    /*
    } else if ( particleName == "carbon" ) {
      pPhysicsHelper->RegisterProcess(new G4DNAIonisation("carbon_G4DNAIonisation"), particle);

    } else if ( particleName == "nitrogen" ) {
      pPhysicsHelper->RegisterProcess(new G4DNAIonisation("nitrogen_G4DNAIonisation"), particle);

    } else if ( particleName == "oxygen" ) {
      pPhysicsHelper->RegisterProcess(new G4DNAIonisation("oxygen_G4DNAIonisation"), particle);

    } else if ( particleName == "iron" ) {
      pPhysicsHelper->RegisterProcess(new G4DNAIonisation("iron_G4DNAIonisation"), particle);
    */ 

    }
    
    // Warning : the following particles and processes are needed by EM Physics builders
    // They are taken from the default Livermore Physics list
    // These particles are currently not handled by Geant4-DNA
    
      // e+
      
    else if (particleName == "e+") {
      
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);      

      pPhysicsHelper->RegisterProcess(msc, pParticle);
      pPhysicsHelper->RegisterProcess(eIoni, pParticle);
      pPhysicsHelper->RegisterProcess(new G4eBremsstrahlung(), pParticle);
      pPhysicsHelper->RegisterProcess(new G4eplusAnnihilation(), pParticle);

    } else if (particleName == "gamma") {
    
      // photoelectric effect - Livermore model only
      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      thePhotoElectricEffect->SetEmModel(new G4LivermorePhotoElectricModel());
      pPhysicsHelper->RegisterProcess(thePhotoElectricEffect, pParticle);

      // Compton scattering - Livermore model only
      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      theComptonScattering->SetEmModel(new G4LivermoreComptonModel());
      pPhysicsHelper->RegisterProcess(theComptonScattering, pParticle);

      // gamma conversion - Livermore model below 80 GeV
      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      theGammaConversion->SetEmModel(new G4LivermoreGammaConversionModel());
      pPhysicsHelper->RegisterProcess(theGammaConversion, pParticle);

      // default Rayleigh scattering is Livermore
      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      pPhysicsHelper->RegisterProcess(theRayleigh, pParticle);
    }    
    // Warning : end of particles and processes are needed by EM Physics builders 
  }

  // Deexcitation
  //
  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
