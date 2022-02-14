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
/// \file ElectromagneticPhysics.cc
/// \brief Implementation of the ElectromagneticPhysics class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ElectromagneticPhysics.hh"

#include "G4BuilderType.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4PhysicsListHelper.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4RayleighScattering.hh"

#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hMultipleScattering.hh"
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

ElectromagneticPhysics::ElectromagneticPhysics(const G4String& name)
   :  G4VPhysicsConstructor(name)
{
    SetPhysicsType(bElectromagnetic);

    G4EmParameters* param = G4EmParameters::Instance();
    param->SetDefaults();
    param->SetStepFunction(0.2, 100*um);
    param->SetStepFunctionMuHad(0.1, 10*um);
    param->SetStepFunctionLightIons(0.1, 10*um);
    param->SetStepFunctionIons(0.1, 1*um);
    param->SetDeexcitationIgnoreCut(true);      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ElectromagneticPhysics::~ElectromagneticPhysics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ElectromagneticPhysics::ConstructProcess()
{
  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  
  // Add standard EM Processes
  //
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    
      ph->RegisterProcess(new G4RayleighScattering,  particle);      
      ph->RegisterProcess(new G4PhotoElectricEffect, particle);
      ph->RegisterProcess(new G4ComptonScattering,   particle);
      ph->RegisterProcess(new G4GammaConversion,     particle);
       
    } else if (particleName == "e-") {
    
      ph->RegisterProcess(new G4eMultipleScattering(), particle);
      ph->RegisterProcess(new G4eIonisation,           particle);
      ph->RegisterProcess(new G4eBremsstrahlung(),     particle);
            
    } else if (particleName == "e+") {
    
      ph->RegisterProcess(new G4eMultipleScattering(), particle);
      ph->RegisterProcess(new G4eIonisation,           particle);
      ph->RegisterProcess(new G4eBremsstrahlung(),     particle);
      ph->RegisterProcess(new G4eplusAnnihilation(),   particle);
                  
    } else if (particleName == "mu+" || 
               particleName == "mu-"    ) {

      ph->RegisterProcess(new G4MuMultipleScattering(), particle);
      ph->RegisterProcess(new G4MuIonisation,           particle);
      ph->RegisterProcess(new G4MuBremsstrahlung(),     particle);
      ph->RegisterProcess(new G4MuPairProduction(),     particle);
                   
    } else if( particleName == "proton" ||
               particleName == "pi-" ||
               particleName == "pi+"    ) {

      ph->RegisterProcess(new G4hMultipleScattering(), particle);      
      ph->RegisterProcess(new G4hIonisation,           particle);
     
    } else if( particleName == "alpha" || 
               particleName == "He3"    ) {

      ph->RegisterProcess(new G4hMultipleScattering(), particle);
      ph->RegisterProcess(new G4ionIonisation,         particle);
      ph->RegisterProcess(new G4NuclearStopping(),     particle);
            
    } else if( particleName == "GenericIon" ) {

      ph->RegisterProcess(new G4hMultipleScattering(), particle);
      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ph->RegisterProcess(ionIoni,                    particle);
      ph->RegisterProcess(new G4NuclearStopping(),    particle);
      
    } else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) && 
               (particle->GetParticleName() != "chargedgeantino")) {
               
      //all others charged particles except geantino
      ph->RegisterProcess(new G4hMultipleScattering(), particle);
      ph->RegisterProcess(new G4hIonisation(),         particle);
    }
  }
 
  // Deexcitation
  //
  G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
  G4LossTableManager::Instance()->SetAtomDeexcitation(de);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

