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

#include "G4EmDNAPhysics_stationary_option2.hh"

#include "G4SystemOfUnits.hh"

#include "G4DNAGenericIonsManager.hh"

// *** Processes and models for Geant4-DNA

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
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_stationary_option2);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysics_stationary_option2::G4EmDNAPhysics_stationary_option2(G4int ver)
  : G4VPhysicsConstructor("G4EmDNAPhysics_stationary_option2"), verbose(ver)
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

G4EmDNAPhysics_stationary_option2::G4EmDNAPhysics_stationary_option2(G4int ver, 
const G4String&)
: G4VPhysicsConstructor("G4EmDNAPhysics_stationary_option2"), verbose(ver)
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

G4EmDNAPhysics_stationary_option2::~G4EmDNAPhysics_stationary_option2()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysics_stationary_option2::ConstructParticle()
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

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysics_stationary_option2::ConstructProcess()
{
  if(verbose > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() )
  {
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    if (particleName == "e-") {

      // *** Elastic scattering (two alternative models available) ***
      
      G4DNAElastic* theDNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
      theDNAElasticProcess->SetEmModel(new G4DNAChampionElasticModel());
      
      // or alternative model
      //theDNAElasticProcess
      //->SetEmModel(new G4DNAScreenedRutherfordElasticModel());
      
      ph->RegisterProcess(theDNAElasticProcess, particle);

      // *** Excitation ***

      G4DNAExcitation* theDNAExcitationProcess = 
       new G4DNAExcitation("e-_G4DNAExcitation");
      theDNAExcitationProcess->SetEmModel(new G4DNABornExcitationModel());
      ((G4DNABornExcitationModel*)(theDNAExcitationProcess->EmModel()))
       ->SelectStationary(true);
      ph->RegisterProcess(theDNAExcitationProcess, particle);

      // *** Ionisation ***

      G4DNAIonisation* theDNAIonisationProcess = 
       new G4DNAIonisation("e-_G4DNAIonisation");
      theDNAIonisationProcess->SetEmModel(new G4DNABornIonisationModel());
      ((G4DNABornIonisationModel*)(theDNAIonisationProcess->EmModel()))
       ->SelectStationary(true);
      //
      ((G4DNABornIonisationModel*)(theDNAIonisationProcess->EmModel()))
       ->SelectFasterComputation(true);
      //
      ph->RegisterProcess(theDNAIonisationProcess, particle);

      // *** Vibrational excitation ***
 
      G4DNAVibExcitation* theDNAVibExcitationProcess = 
       new G4DNAVibExcitation("e-_G4DNAVibExcitation");
      theDNAVibExcitationProcess->SetEmModel(new G4DNASancheExcitationModel());
      ((G4DNASancheExcitationModel*)(theDNAVibExcitationProcess->EmModel()))
       ->SelectStationary(true);
      ph->RegisterProcess(theDNAVibExcitationProcess, particle);
          
      // *** Attachment ***

      G4DNAAttachment* theDNAAttachmentProcess = 
       new G4DNAAttachment("e-_G4DNAAttachment");
      theDNAAttachmentProcess->SetEmModel(new G4DNAMeltonAttachmentModel());
      ((G4DNAMeltonAttachmentModel*)(theDNAAttachmentProcess->EmModel()))
       ->SelectStationary(true);
      ph->RegisterProcess(theDNAAttachmentProcess, particle);
   
    } else if ( particleName == "proton" ) {
      
      // *** Elastic ***
      
      G4DNAElastic* theDNAElasticProcess = 
       new G4DNAElastic("proton_G4DNAElastic");
      theDNAElasticProcess->SetEmModel(new G4DNAIonElasticModel());
      ((G4DNAIonElasticModel*)(theDNAElasticProcess->EmModel()))
       ->SelectStationary(true);
      ph->RegisterProcess(theDNAElasticProcess, particle);
      
      // *** Excitation ***
      
      G4DNAExcitation* theDNAExcitationProcess = 
       new G4DNAExcitation("proton_G4DNAExcitation");
     
      theDNAExcitationProcess->SetEmModel
       (new G4DNAMillerGreenExcitationModel(),1);
      theDNAExcitationProcess->SetEmModel
       (new G4DNABornExcitationModel(),2);
      
      ((G4DNAMillerGreenExcitationModel*)(theDNAExcitationProcess->EmModel(1)))
       ->SetLowEnergyLimit(10*eV);
      ((G4DNAMillerGreenExcitationModel*)(theDNAExcitationProcess->EmModel(1)))
       ->SetHighEnergyLimit(500*keV);
      ((G4DNAMillerGreenExcitationModel*)(theDNAExcitationProcess->EmModel(1)))
       ->SelectStationary(true);
      
      ((G4DNABornExcitationModel*)(theDNAExcitationProcess->EmModel(2)))
       ->SetLowEnergyLimit(500*keV);
      ((G4DNABornExcitationModel*)(theDNAExcitationProcess->EmModel(2)))
       ->SetHighEnergyLimit(100*MeV);
      ((G4DNABornExcitationModel*)(theDNAExcitationProcess->EmModel(2)))
       ->SelectStationary(true);
      
      ph->RegisterProcess(theDNAExcitationProcess, particle);

      // *** Ionisation ***
      
      G4DNAIonisation* theDNAIonisationProcess = 
       new G4DNAIonisation("proton_G4DNAIonisation");

      theDNAIonisationProcess->SetEmModel(
       new G4DNARuddIonisationExtendedModel,1);
      theDNAIonisationProcess->SetEmModel(
       new G4DNABornIonisationModel,2);
      
      ((G4DNARuddIonisationExtendedModel*)
       (theDNAIonisationProcess->EmModel(1)))->SetLowEnergyLimit(0*eV);
      ((G4DNARuddIonisationExtendedModel*)
       (theDNAIonisationProcess->EmModel(1)))->SetHighEnergyLimit(500*keV);
      ((G4DNARuddIonisationExtendedModel*)
       (theDNAIonisationProcess->EmModel(1)))->SelectStationary(true);
      
      ((G4DNABornIonisationModel*)
       (theDNAIonisationProcess->EmModel(2)))->SetLowEnergyLimit(500*keV);
      ((G4DNABornIonisationModel*)
       (theDNAIonisationProcess->EmModel(2)))->SetHighEnergyLimit(100*MeV);
      ((G4DNABornIonisationModel*)
       (theDNAIonisationProcess->EmModel(2)))->SelectStationary(true);
      //
      ((G4DNABornIonisationModel*)
       (theDNAIonisationProcess->EmModel(2)))->SelectFasterComputation(true);
      //

      ph->RegisterProcess(theDNAIonisationProcess, particle);

      // *** Charge decrease ***

      G4DNAChargeDecrease* theDNAChargeDecreaseProcess = 
       new G4DNAChargeDecrease("proton_G4DNAChargeDecrease");
      theDNAChargeDecreaseProcess->SetEmModel(
       new G4DNADingfelderChargeDecreaseModel());
      ((G4DNADingfelderChargeDecreaseModel*)
       (theDNAChargeDecreaseProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAChargeDecreaseProcess, particle);

    } else if ( particleName == "hydrogen" ) {
      
      // *** Elastic ***
      
      G4DNAElastic* theDNAElasticProcess = 
       new G4DNAElastic("hydrogen_G4DNAElastic");
      theDNAElasticProcess->SetEmModel(new G4DNAIonElasticModel());
      ((G4DNAIonElasticModel*)(theDNAElasticProcess->EmModel()))
       ->SelectStationary(true);
      ph->RegisterProcess(theDNAElasticProcess, particle);
      
      // *** Excitation ***

      G4DNAExcitation* theDNAExcitationProcess = 
       new G4DNAExcitation("hydrogen_G4DNAExcitation");
      theDNAExcitationProcess->SetEmModel(
       new G4DNAMillerGreenExcitationModel());
      ((G4DNAMillerGreenExcitationModel*)(theDNAExcitationProcess->EmModel()))
       ->SelectStationary(true);
      ph->RegisterProcess(theDNAExcitationProcess, particle);
      
      // *** Ionisation ***

      G4DNAIonisation* theDNAIonisationProcess = 
       new G4DNAIonisation("hydrogen_G4DNAIonisation");
      theDNAIonisationProcess->SetEmModel(
       new G4DNARuddIonisationExtendedModel());
      ((G4DNARuddIonisationExtendedModel*)
       (theDNAIonisationProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAIonisationProcess, particle);
      
      // *** Charge increase ***

      G4DNAChargeIncrease* theDNAChargeIncreaseProcess = 
       new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease");
      theDNAChargeIncreaseProcess->SetEmModel(
       new G4DNADingfelderChargeIncreaseModel());
      ((G4DNADingfelderChargeIncreaseModel*)
       (theDNAChargeIncreaseProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAChargeIncreaseProcess, particle);

    } else if ( particleName == "alpha" ) {
      
      // *** Elastic ***
      
      G4DNAElastic* theDNAElasticProcess = 
       new G4DNAElastic("alpha_G4DNAElastic");
      theDNAElasticProcess->SetEmModel(new G4DNAIonElasticModel());
      ((G4DNAIonElasticModel*)
       (theDNAElasticProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAElasticProcess, particle);
      
      // *** Excitation ***

      G4DNAExcitation* theDNAExcitationProcess = 
       new G4DNAExcitation("alpha_G4DNAExcitation");
      theDNAExcitationProcess->SetEmModel(
       new G4DNAMillerGreenExcitationModel());
      ((G4DNAMillerGreenExcitationModel*)
       (theDNAExcitationProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAExcitationProcess, particle);
            
      // *** Ionisation ***

      G4DNAIonisation* theDNAIonisationProcess = 
       new G4DNAIonisation("alpha_G4DNAIonisation");
      theDNAIonisationProcess->SetEmModel(
       new G4DNARuddIonisationExtendedModel());
      ((G4DNARuddIonisationExtendedModel*)
       (theDNAIonisationProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAIonisationProcess, particle);

      // *** Charge decrease ***

      G4DNAChargeDecrease* theDNAChargeDecreaseProcess = 
       new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease");
      theDNAChargeDecreaseProcess->SetEmModel(
       new G4DNADingfelderChargeDecreaseModel());
      ((G4DNADingfelderChargeDecreaseModel*)
       (theDNAChargeDecreaseProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAChargeDecreaseProcess, particle);

    } else if ( particleName == "alpha+" ) {
      
      // *** Elastic ***
      
      G4DNAElastic* theDNAElasticProcess = 
       new G4DNAElastic("alpha+_G4DNAElastic");
      theDNAElasticProcess->SetEmModel(new G4DNAIonElasticModel());
      ((G4DNAIonElasticModel*)
       (theDNAElasticProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAElasticProcess, particle);
      
      // *** Excitation ***

      G4DNAExcitation* theDNAExcitationProcess = 
       new G4DNAExcitation("alpha+_G4DNAExcitation");
      theDNAExcitationProcess->SetEmModel(
       new G4DNAMillerGreenExcitationModel());
      ((G4DNAMillerGreenExcitationModel*)
       (theDNAExcitationProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAExcitationProcess, particle);
            
      // *** Ionisation ***

      G4DNAIonisation* theDNAIonisationProcess = 
       new G4DNAIonisation("alpha+_G4DNAIonisation");
      theDNAIonisationProcess->SetEmModel(
       new G4DNARuddIonisationExtendedModel());
      ((G4DNARuddIonisationExtendedModel*)
       (theDNAIonisationProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAIonisationProcess, particle);

      // *** Charge decrease ***

      G4DNAChargeDecrease* theDNAChargeDecreaseProcess = 
       new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease");
      theDNAChargeDecreaseProcess->SetEmModel(
       new G4DNADingfelderChargeDecreaseModel());
      ((G4DNADingfelderChargeDecreaseModel*)
       (theDNAChargeDecreaseProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAChargeDecreaseProcess, particle);

      // *** Charge increase ***

      G4DNAChargeIncrease* theDNAChargeIncreaseProcess = 
       new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease");
      theDNAChargeIncreaseProcess->SetEmModel(
       new G4DNADingfelderChargeIncreaseModel());
      ((G4DNADingfelderChargeIncreaseModel*)
       (theDNAChargeIncreaseProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAChargeIncreaseProcess, particle);

    } else if ( particleName == "helium" ) {
      
      // *** Elastic ***
      
      G4DNAElastic* theDNAElasticProcess = 
       new G4DNAElastic("helium_G4DNAElastic");
      theDNAElasticProcess->SetEmModel(
       new G4DNAIonElasticModel());
      ((G4DNAIonElasticModel*)
       (theDNAElasticProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAElasticProcess, particle);
      
      // *** Excitation ***

      G4DNAExcitation* theDNAExcitationProcess = 
       new G4DNAExcitation("helium_G4DNAExcitation");
      theDNAExcitationProcess->SetEmModel(
       new G4DNAMillerGreenExcitationModel());
      ((G4DNAMillerGreenExcitationModel*)
       (theDNAExcitationProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAExcitationProcess, particle);
            
      // *** Ionisation ***

      G4DNAIonisation* theDNAIonisationProcess = 
       new G4DNAIonisation("helium_G4DNAIonisation");
      theDNAIonisationProcess->SetEmModel(
       new G4DNARuddIonisationExtendedModel());
      ((G4DNARuddIonisationExtendedModel*)
       (theDNAIonisationProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAIonisationProcess, particle);

      // *** Charge increase ***

      G4DNAChargeIncrease* theDNAChargeIncreaseProcess = 
       new G4DNAChargeIncrease("helium_G4DNAChargeIncrease");
      theDNAChargeIncreaseProcess->SetEmModel(
       new G4DNADingfelderChargeIncreaseModel());
      ((G4DNADingfelderChargeIncreaseModel*)
       (theDNAChargeIncreaseProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAChargeIncreaseProcess, particle);
    
    } else if ( particleName == "GenericIon" ) {

      // *** Ionisation ***

      G4DNAIonisation* theDNAIonisationProcess = 
       new G4DNAIonisation("GenericIon_G4DNAIonisation");
      theDNAIonisationProcess->SetEmModel(
       new G4DNARuddIonisationExtendedModel());
      ((G4DNARuddIonisationExtendedModel*)
       (theDNAIonisationProcess->EmModel()))->SelectStationary(true);
      ph->RegisterProcess(theDNAIonisationProcess, particle);

    }
    
    // Warning : the following particles and processes are needed by EM Physics 
    // builders
    // They are taken from the default Livermore Physics list
    // These particles are currently not handled by Geant4-DNA
    
      // e+
      
    else if (particleName == "e+") {

      // Identical to G4EmStandardPhysics_stationary
      
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);      

      ph->RegisterProcess(msc, particle);
      ph->RegisterProcess(eIoni, particle);
      ph->RegisterProcess(new G4eBremsstrahlung(), particle);
      ph->RegisterProcess(new G4eplusAnnihilation(), particle);

    } else if (particleName == "gamma") {
    
      G4double LivermoreHighEnergyLimit = GeV;

      G4PhotoElectricEffect* thePhotoElectricEffect = 
       new G4PhotoElectricEffect();
      G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel = 
       new G4LivermorePhotoElectricModel();
      theLivermorePhotoElectricModel
       ->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      thePhotoElectricEffect->AddEmModel(0, theLivermorePhotoElectricModel);
      ph->RegisterProcess(thePhotoElectricEffect, particle);

      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      G4LivermoreComptonModel* theLivermoreComptonModel = 
       new G4LivermoreComptonModel();
      theLivermoreComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theComptonScattering->AddEmModel(0, theLivermoreComptonModel);
      ph->RegisterProcess(theComptonScattering, particle);

      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      G4LivermoreGammaConversionModel* theLivermoreGammaConversionModel = 
       new G4LivermoreGammaConversionModel();
      theLivermoreGammaConversionModel
       ->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theGammaConversion->AddEmModel(0, theLivermoreGammaConversionModel);
      ph->RegisterProcess(theGammaConversion, particle);

      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      G4LivermoreRayleighModel* theRayleighModel = 
       new G4LivermoreRayleighModel();
      theRayleighModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theRayleigh->AddEmModel(0, theRayleighModel);
      ph->RegisterProcess(theRayleigh, particle);
    }
    
   // Warning : end of particles and processes are needed by EM Physics build. 
    
  }

  // Deexcitation
  //
  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);
  de->SetFluo(true);
}

