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
// $Id: G4EmDNAPhysics.cc,v 1.8 2010-11-25 07:44:55 sincerti Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "G4EmDNAPhysics.hh"

#include "G4ProcessManager.hh"
#include "G4DNAGenericIonsManager.hh"

// *** Processes and models for Geant4-DNA

#include "G4DNAElastic.hh"
#include "G4DNAChampionElasticModel.hh"
#include "G4DNAScreenedRutherfordElasticModel.hh"

#include "G4DNAExcitation.hh"
#include "G4DNAAttachment.hh"
#include "G4DNAVibExcitation.hh"
#include "G4DNAIonisation.hh"
#include "G4DNAChargeDecrease.hh"
#include "G4DNAChargeIncrease.hh"

// particles

#include "G4Electron.hh"
#include "G4Proton.hh"

// Warning : the following is needed in order to use EM Physics builders
// e+
#include "G4Positron.hh"
#include "G4eMultipleScattering.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
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
// end of warning

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysics::G4EmDNAPhysics(G4int ver)
  : G4VPhysicsConstructor("G4EmDNAPhysics"), verbose(ver)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysics::G4EmDNAPhysics(G4int ver, const G4String&)
  : G4VPhysicsConstructor("G4EmDNAPhysics"), verbose(ver)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDNAPhysics::~G4EmDNAPhysics()
{}

#include "G4GenericIon.hh"

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

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDNAPhysics::ConstructProcess()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "e-") {

      // *** Elastic scattering (two alternative models available) ***
      
      G4DNAElastic* theDNAElasticProcess = new G4DNAElastic();
      // theDNAElasticProcess->SetModel(new G4DNAChampionElasticModel());
      
      // or alternative model
      theDNAElasticProcess->SetModel(new G4DNAScreenedRutherfordElasticModel());
      
      pmanager->AddDiscreteProcess(theDNAElasticProcess);

      // *** Excitation ***
      pmanager->AddDiscreteProcess(new G4DNAExcitation());

      // *** Ionisation ***
      pmanager->AddDiscreteProcess(new G4DNAIonisation());

      // *** Vibrational excitation ***
      pmanager->AddDiscreteProcess(new G4DNAVibExcitation());
      
      // *** Attachment ***
      pmanager->AddDiscreteProcess(new G4DNAAttachment()); 

    
    } else if ( particleName == "proton" ) {
      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeDecrease());

    } else if ( particleName == "hydrogen" ) {
      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeIncrease());

    } else if ( particleName == "alpha" ) {
      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeDecrease());

    } else if ( particleName == "alpha+" ) {
      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeDecrease());
      pmanager->AddDiscreteProcess(new G4DNAChargeIncrease());

    } else if ( particleName == "helium" ) {
      pmanager->AddDiscreteProcess(new G4DNAExcitation());
      pmanager->AddDiscreteProcess(new G4DNAIonisation());
      pmanager->AddDiscreteProcess(new G4DNAChargeIncrease());

    }
    
    // Warning : the following particles and processes are needed by EM Physics builders
    // They are taken from the default Livermore Physics list
    // These particles are currently not handled by Geant4-DNA
    
      // e+
      
      else if (particleName == "e+") {

      // Identical to G4EmStandardPhysics_option3
      
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->SetStepLimitType(fUseDistanceToBoundary);
      msc->AddEmModel(0, new G4GoudsmitSaundersonMscModel());
      pmanager->AddProcess(msc,                   -1, 1, 1);

      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);      

      pmanager->AddProcess(eIoni,                 -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung, -1,-3, 3);      
      pmanager->AddProcess(new G4eplusAnnihilation,0,-1, 4);

    }  else if (particleName == "gamma") {
    
      G4double LivermoreHighEnergyLimit = GeV;

      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel = 
	new G4LivermorePhotoElectricModel();
      theLivermorePhotoElectricModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      thePhotoElectricEffect->AddEmModel(0, theLivermorePhotoElectricModel);
      pmanager->AddDiscreteProcess(thePhotoElectricEffect);

      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      G4LivermoreComptonModel* theLivermoreComptonModel = 
	new G4LivermoreComptonModel();
      theLivermoreComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theComptonScattering->AddEmModel(0, theLivermoreComptonModel);
      pmanager->AddDiscreteProcess(theComptonScattering);

      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      G4LivermoreGammaConversionModel* theLivermoreGammaConversionModel = 
	new G4LivermoreGammaConversionModel();
      theLivermoreGammaConversionModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theGammaConversion->AddEmModel(0, theLivermoreGammaConversionModel);
      pmanager->AddDiscreteProcess(theGammaConversion);

      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      G4LivermoreRayleighModel* theRayleighModel = new G4LivermoreRayleighModel();
      theRayleighModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theRayleigh->AddEmModel(0, theRayleighModel);
      pmanager->AddDiscreteProcess(theRayleigh);
    }
    
   // Warning : end of particles and processes are needed by EM Physics builders 
    
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
