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
/// \file electromagnetic/TestEm18/src/PhysListEmStandard.cc
/// \brief Implementation of the PhysListEmStandard class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysListEmStandard.hh"

#include "G4BuilderType.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4PhysicsListHelper.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"
#include "G4KleinNishinaModel.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"

#include "G4ionIonisation.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"

#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::PhysListEmStandard(const G4String& name)
   :  G4VPhysicsConstructor(name)
{
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetDefaults();
  param->SetMinEnergy(10*eV);
  param->SetMaxEnergy(10*TeV);
  param->SetNumberOfBinsPerDecade(10);
  param->SetBuildCSDARange(true);
  param->SetMaxEnergyForCSDARange(10*TeV);
  SetPhysicsType(bElectromagnetic);
  
  param->SetVerbose(0);
  param->Dump();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard::~PhysListEmStandard()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandard::ConstructProcess()
{
  G4PhysicsListHelper* list = G4PhysicsListHelper::GetPhysicsListHelper();
  
  // Add standard EM Processes
  //
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {

      ////list->RegisterProcess(new G4RayleighScattering, particle);
      list->RegisterProcess(new G4PhotoElectricEffect, particle);
      G4ComptonScattering* compt   = new G4ComptonScattering;
      compt->SetEmModel(new G4KleinNishinaModel());
      list->RegisterProcess(compt, particle);
      list->RegisterProcess(new G4GammaConversion, particle);
     
    } else if (particleName == "e-") {

      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.1, 100*um);      
      list->RegisterProcess(eIoni,                   particle);
      list->RegisterProcess(new G4eBremsstrahlung(), particle);
            
    } else if (particleName == "e+") {

      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.1, 100*um);      
      list->RegisterProcess(eIoni,                     particle);
      list->RegisterProcess(new G4eBremsstrahlung(),   particle);
      list->RegisterProcess(new G4eplusAnnihilation(), particle);
                  
    } else if (particleName == "mu+" || 
               particleName == "mu-"    ) {

      G4MuIonisation* muIoni = new G4MuIonisation();
      muIoni->SetStepFunction(0.1, 50*um);      
      list->RegisterProcess(muIoni,                   particle);
      list->RegisterProcess(new G4MuBremsstrahlung(), particle);
      list->RegisterProcess(new G4MuPairProduction(), particle);
                   
    } else if( particleName == "proton" ||
               particleName == "pi-" ||
               particleName == "pi+"    ) {
  
      G4hIonisation* hIoni = new G4hIonisation();
      hIoni->SetStepFunction(0.1, 20*um);
      list->RegisterProcess(hIoni,                   particle);
      list->RegisterProcess(new G4hBremsstrahlung(), particle);
      list->RegisterProcess(new G4hPairProduction(), particle);
     
    } else if( particleName == "alpha" || 
               particleName == "He3"    ) {

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetStepFunction(0.1, 1*um);
      list->RegisterProcess(ionIoni,                 particle);
      list->RegisterProcess(new G4NuclearStopping(), particle);
            
    } else if( particleName == "GenericIon" ) {
 
      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ionIoni->SetStepFunction(0.1, 1*um);
      list->RegisterProcess(ionIoni,                 particle);
      list->RegisterProcess(new G4NuclearStopping(), particle);
      
    } else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) && 
               (particle->GetParticleName() != "chargedgeantino")) {
               
      //all others charged particles except geantino
      list->RegisterProcess(new G4hIonisation(),    particle);
    }
  }
    
  // Deexcitation
  //
  G4VAtomDeexcitation* deex = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(deex);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

