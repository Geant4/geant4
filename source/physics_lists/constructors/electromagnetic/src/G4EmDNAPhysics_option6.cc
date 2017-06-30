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
// Based on the work of M. Terrissol and M. C. Bordage
//
// Users are requested to cite the following papers:
// - M. Terrissol, A. Baudre, Radiat. Prot. Dosim. 31 (1990) 175-177
// - M.C. Bordage, J. Bordes, S. Edel, M. Terrissol, X. Franceries, 
//   M. Bardies, N. Lampe, S. Incerti, Phys. Med. 32 (2016) 1833-1840
//
// Authors of this class: 
// M.C. Bordage, M. Terrissol, S. Edel, J. Bordes, S. Incerti
//
// 15.01.2014: creation
//

#include "G4EmDNAPhysics_option6.hh"

#include "G4SystemOfUnits.hh"

#include "G4DNAGenericIonsManager.hh"

// *** Processes and models for Geant4-DNA

#include "G4DNAElectronSolvation.hh"
#include "G4DNAOneStepThermalizationModel.hh"

#include "G4DNAElastic.hh"
#include "G4DNACPA100ElasticModel.hh"

#include "G4DNAIonisation.hh"
#include "G4DNACPA100IonisationModel.hh"

#include "G4DNAExcitation.hh"
#include "G4DNACPA100ExcitationModel.hh"

#include "G4DNAAttachment.hh"
#include "G4DNAVibExcitation.hh"

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
G4_DECLARE_PHYSCONSTR_FACTORY(G4EmDNAPhysics_option6);


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysics_option6::G4EmDNAPhysics_option6(G4int ver, const G4String&)
  : G4VPhysicsConstructor("G4EmDNAPhysics_option6"), verbose(ver)
{
//  G4EmParameters::Instance()->SetDefaults();
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetFluo(true);  
  param->SetAuger(true);  
  param->SetAugerCascade(true);  
  param->SetDeexcitationIgnoreCut(true);

  SetPhysicsType(bElectromagnetic);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysics_option6::~G4EmDNAPhysics_option6()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysics_option6::ConstructParticle()
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

void G4EmDNAPhysics_option6::ConstructProcess()
{
  if(verbose > 1) {
    G4cout << "### " << GetPhysicsName() << " Construct Processes " << G4endl;
  }
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();

  //OLD
  //aParticleIterator->reset();
  //while( (*aParticleIterator)() )

  auto myParticleIterator=GetParticleIterator();
  myParticleIterator->reset();
  while( (*myParticleIterator)() )

  {
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4String particleName = particle->GetParticleName();

    if (particleName == "e-") {

      // *** Solvation ***
      G4DNAElectronSolvation* solvation = 
        new G4DNAElectronSolvation("e-_G4DNAElectronSolvation");
      G4DNAOneStepThermalizationModel* therm =
       new G4DNAOneStepThermalizationModel();
      therm->SetHighEnergyLimit(11.*eV); // limit of the CPA100 elastic model
      solvation->SetEmModel(therm);
      ph->RegisterProcess(solvation, particle);

      // *** Elastic scattering (two alternative models available) ***      
      G4DNAElastic* theDNAElasticProcess = new G4DNAElastic("e-_G4DNAElastic");
      theDNAElasticProcess->SetEmModel(new G4DNACPA100ElasticModel());
      ph->RegisterProcess(theDNAElasticProcess, particle);

      // *** Excitation ***
      G4DNAExcitation* theDNAExcitationProcess = new G4DNAExcitation("e-_G4DNAExcitation");
      theDNAExcitationProcess->SetEmModel(new G4DNACPA100ExcitationModel());
      ph->RegisterProcess(theDNAExcitationProcess, particle);

      // *** Ionisation ***
      G4DNAIonisation* theDNAIonisationProcess = new G4DNAIonisation("e-_G4DNAIonisation");
      theDNAIonisationProcess->SetEmModel(new G4DNACPA100IonisationModel());
      ph->RegisterProcess(theDNAIonisationProcess, particle);

      // *** Vibrational excitation ***
      //ph->RegisterProcess(new G4DNAVibExcitation("e-_G4DNAVibExcitation"), particle);
      
      // *** Attachment ***
      //ph->RegisterProcess(new G4DNAAttachment("e-_G4DNAAttachment"), particle); 
    
    } else if ( particleName == "proton" ) {
      ph->RegisterProcess(new G4DNAElastic("proton_G4DNAElastic"), particle);
      ph->RegisterProcess(new G4DNAExcitation("proton_G4DNAExcitation"), particle);
      ph->RegisterProcess(new G4DNAIonisation("proton_G4DNAIonisation"), particle);
      ph->RegisterProcess(new G4DNAChargeDecrease("proton_G4DNAChargeDecrease"), particle);

    } else if ( particleName == "hydrogen" ) {
      ph->RegisterProcess(new G4DNAElastic("hydrogen_G4DNAElastic"), particle);
      ph->RegisterProcess(new G4DNAExcitation("hydrogen_G4DNAExcitation"), particle);
      ph->RegisterProcess(new G4DNAIonisation("hydrogen_G4DNAIonisation"), particle);
      ph->RegisterProcess(new G4DNAChargeIncrease("hydrogen_G4DNAChargeIncrease"), particle);

    } else if ( particleName == "alpha" ) {
      ph->RegisterProcess(new G4DNAElastic("alpha_G4DNAElastic"), particle);
      ph->RegisterProcess(new G4DNAExcitation("alpha_G4DNAExcitation"), particle);
      ph->RegisterProcess(new G4DNAIonisation("alpha_G4DNAIonisation"), particle);
      ph->RegisterProcess(new G4DNAChargeDecrease("alpha_G4DNAChargeDecrease"), particle);

    } else if ( particleName == "alpha+" ) {
      ph->RegisterProcess(new G4DNAElastic("alpha+_G4DNAElastic"), particle);
      ph->RegisterProcess(new G4DNAExcitation("alpha+_G4DNAExcitation"), particle);
      ph->RegisterProcess(new G4DNAIonisation("alpha+_G4DNAIonisation"), particle);
      ph->RegisterProcess(new G4DNAChargeDecrease("alpha+_G4DNAChargeDecrease"), particle);
      ph->RegisterProcess(new G4DNAChargeIncrease("alpha+_G4DNAChargeIncrease"), particle);

    } else if ( particleName == "helium" ) {
      ph->RegisterProcess(new G4DNAElastic("helium_G4DNAElastic"), particle);
      ph->RegisterProcess(new G4DNAExcitation("helium_G4DNAExcitation"), particle);
      ph->RegisterProcess(new G4DNAIonisation("helium_G4DNAIonisation"), particle);
      ph->RegisterProcess(new G4DNAChargeIncrease("helium_G4DNAChargeIncrease"), particle);
    
    } else if ( particleName == "GenericIon" ) {
      ph->RegisterProcess(new G4DNAIonisation("GenericIon_G4DNAIonisation"), particle);
    }
    
    // Warning : the following particles and processes are needed by EM Physics builders
    // They are taken from the default Livermore Physics list
    // These particles are currently not handled by Geant4-DNA
    
      // e+
      
    else if (particleName == "e+") {

      // Identical to G4EmStandardPhysics_option3
      
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

      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel = 
	new G4LivermorePhotoElectricModel();
      theLivermorePhotoElectricModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
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
      theLivermoreGammaConversionModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theGammaConversion->AddEmModel(0, theLivermoreGammaConversionModel);
      ph->RegisterProcess(theGammaConversion, particle);

      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      G4LivermoreRayleighModel* theRayleighModel = new G4LivermoreRayleighModel();
      theRayleighModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theRayleigh->AddEmModel(0, theRayleighModel);
      ph->RegisterProcess(theRayleigh, particle);

    }
    
   // Warning : end of particles and processes are needed by EM Physics builders 
    
  }

  // Deexcitation
  //

  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
