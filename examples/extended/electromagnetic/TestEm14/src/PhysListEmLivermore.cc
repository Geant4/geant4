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
/// \file PhysListEmLivermore.cc
/// \brief Implementation of the PhysListEmLivermore class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysListEmLivermore.hh"
#include "G4BuilderType.hh"
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

#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"
#include "G4UniversalFluctuation.hh"

#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"

// e+

#include "G4eplusAnnihilation.hh"

// mu

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

// hadrons, ions

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

// deexcitation

#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmLivermore::PhysListEmLivermore(const G4String& name)
  :  G4VPhysicsConstructor(name)
{
    G4EmParameters* param = G4EmParameters::Instance();
    param->SetDefaults();
    param->SetMinEnergy(10*eV);
    param->SetMaxEnergy(10*TeV);
    param->SetNumberOfBinsPerDecade(10);
    param->SetBuildCSDARange(true);
  
    param->SetVerbose(0);
    param->Dump();
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmLivermore::ConstructProcess()
{
  // Add Livermore EM Processes
  
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();    
    G4String particleName = particle->GetParticleName();

    //Applicability range for Livermore models
    //for higher energies, the Standard models are used   
    G4double highEnergyLimit = 1*GeV;
         
    if (particleName == "gamma") {
      // gamma         

      G4PhotoElectricEffect* phot = new G4PhotoElectricEffect();
      G4LivermorePhotoElectricModel* 
      photModel = new G4LivermorePhotoElectricModel();
      photModel->SetHighEnergyLimit(highEnergyLimit);
      phot->AddEmModel(0, photModel);
      ///list->RegisterProcess(phot, particle);
      pmanager->AddDiscreteProcess(phot);
      
      G4ComptonScattering* compt = new G4ComptonScattering();
      G4LivermoreComptonModel* 
      comptModel = new G4LivermoreComptonModel();
      comptModel->SetHighEnergyLimit(highEnergyLimit);
      compt->AddEmModel(0, comptModel);
      pmanager->AddDiscreteProcess(compt);

      G4GammaConversion* conv = new G4GammaConversion();
      G4LivermoreGammaConversionModel* 
      convModel = new G4LivermoreGammaConversionModel();
      convModel->SetHighEnergyLimit(highEnergyLimit);
      conv->AddEmModel(0, convModel);
      pmanager->AddDiscreteProcess(conv);
      
      G4RayleighScattering* rayl = new G4RayleighScattering();
      G4LivermoreRayleighModel* 
      raylModel = new G4LivermoreRayleighModel();
      raylModel->SetHighEnergyLimit(highEnergyLimit);
      rayl->AddEmModel(0, raylModel);
      pmanager->AddDiscreteProcess(rayl);
      
    } else if (particleName == "e-") {
      //electron

      G4eIonisation* eIoni = new G4eIonisation();
      G4LivermoreIonisationModel* 
      eIoniModel = new G4LivermoreIonisationModel();
      eIoniModel->SetHighEnergyLimit(highEnergyLimit); 
      eIoni->AddEmModel(0, eIoniModel, new G4UniversalFluctuation() );
      pmanager->AddProcess(eIoni,                   -1,-1, 1);
      
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
      G4LivermoreBremsstrahlungModel* 
      eBremModel = new G4LivermoreBremsstrahlungModel();
      eBremModel->SetHighEnergyLimit(highEnergyLimit);
      eBrem->AddEmModel(0, eBremModel);
      pmanager->AddProcess(eBrem,                   -1,-1, 2);

    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4eIonisation,       -1,-1, 1);
      pmanager->AddProcess(new G4eBremsstrahlung,   -1,-1, 2);
      pmanager->AddProcess(new G4eplusAnnihilation,  0,-1, 3);
      
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon  
      pmanager->AddProcess(new G4MuIonisation,      -1,-1, 1);
      pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1, 2);
      pmanager->AddProcess(new G4MuPairProduction,  -1,-1, 3);       
     
    } else if( particleName == "alpha" || particleName == "GenericIon" ) { 
      pmanager->AddProcess(new G4ionIonisation,     -1,-1, 1);

    } else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) && 
               (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddProcess(new G4hIonisation,       -1,-1, 1);
    }    
  }
  
  // Deexcitation
  //
  G4VAtomDeexcitation* deex = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(deex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

