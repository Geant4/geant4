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
/// \file radioactivedecay/rdecay02/src/exrdmPhysListEmLowEnergy.cc
/// \brief Implementation of the exrdmPhysListEmLowEnergy class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "exrdmPhysListEmLowEnergy.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

// gamma
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"

#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"


// e-
#include "G4eMultipleScattering.hh"
#include "G4UniversalFluctuation.hh"

#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"

// e+
#include "G4eplusAnnihilation.hh"

// mu

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

// hadrons

#include "G4hMultipleScattering.hh"
#include "G4MscStepLimitType.hh"

#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4alphaIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

// msc models
#include "G4UrbanMscModel93.hh"
#include "G4WentzelVIModel.hh"
#include "G4GoudsmitSaundersonMscModel.hh"
#include "G4CoulombScattering.hh"

#include "G4SystemOfUnits.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exrdmPhysListEmLowEnergy::exrdmPhysListEmLowEnergy(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exrdmPhysListEmLowEnergy::~exrdmPhysListEmLowEnergy()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exrdmPhysListEmLowEnergy::ConstructProcess()
{
  // Add EM Processes from G4EmLivermorePhysics builder

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    G4double LivermoreHighEnergyLimit = GeV;
     
    if (particleName == "gamma") {
     
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
      
    } else if (particleName == "e-") {
  
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      //msc->AddEmModel(0, new G4UrbanMscModel93());
      msc->AddEmModel(0, new G4GoudsmitSaundersonMscModel());
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc,                   -1, 1, 1);
      
      // Ionisation
      G4eIonisation* eIoni = new G4eIonisation();
      G4LivermoreIonisationModel* theIoniLivermore = new
        G4LivermoreIonisationModel();
      theIoniLivermore->SetHighEnergyLimit(1*MeV); 
      eIoni->AddEmModel(0, theIoniLivermore, new G4UniversalFluctuation() );
      eIoni->SetStepFunction(0.2, 100*um); //     
      pmanager->AddProcess(eIoni,                 -1, 2, 2);
      
      // Bremsstrahlung
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
      G4LivermoreBremsstrahlungModel* theBremLivermore = new
        G4LivermoreBremsstrahlungModel();
      theBremLivermore->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      eBrem->AddEmModel(0, theBremLivermore);
      pmanager->AddProcess(eBrem, -1,-3, 3);
            
    } else if (particleName == "e+") {

      // Identical to G4EmStandardPhysics_option3
      
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      //msc->AddEmModel(0, new G4UrbanMscModel93());
      msc->AddEmModel(0, new G4GoudsmitSaundersonMscModel());
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc,                   -1, 1, 1);

      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*um);      

      pmanager->AddProcess(eIoni,                 -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung, -1,-3, 3);      
      pmanager->AddProcess(new G4eplusAnnihilation,0,-1, 4);

      
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {

      // Identical to G4EmStandardPhysics_option3
      
      G4MuMultipleScattering* msc = new G4MuMultipleScattering();
      msc->AddEmModel(0, new G4WentzelVIModel());
      pmanager->AddProcess(msc,                       -1, 1, 1);

      G4MuIonisation* muIoni = new G4MuIonisation();
      muIoni->SetStepFunction(0.2, 50*um);          

      pmanager->AddProcess(muIoni,                    -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung,    -1,-3, 3);
      pmanager->AddProcess(new G4MuPairProduction,    -1,-4, 4);
      pmanager->AddDiscreteProcess(new G4CoulombScattering());

    } else if (particleName == "GenericIon") {

      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ionIoni->SetStepFunction(0.1, 10*um);
      pmanager->AddProcess(ionIoni,                   -1, 2, 2);
      pmanager->AddProcess(new G4NuclearStopping(),   -1, 3,-1);

    } else if (particleName == "alpha" ||
               particleName == "He3" ) {

      // Identical to G4EmStandardPhysics_option3
      
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetStepFunction(0.1, 20*um);
      pmanager->AddProcess(ionIoni,                   -1, 2, 2);
      pmanager->AddProcess(new G4NuclearStopping(),   -1, 3,-1);

    } else if (particleName == "pi+" ||
               particleName == "pi-" ||
               particleName == "kaon+" ||
               particleName == "kaon-" ||
               particleName == "proton" ) {

      // Identical to G4EmStandardPhysics_option3
      
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      
      G4hIonisation* hIoni = new G4hIonisation();
      hIoni->SetStepFunction(0.2, 50*um);

      pmanager->AddProcess(hIoni,                     -1, 2, 2);      
      pmanager->AddProcess(new G4hBremsstrahlung,     -1,-3, 3);
      pmanager->AddProcess(new G4hPairProduction,     -1,-4, 4);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

