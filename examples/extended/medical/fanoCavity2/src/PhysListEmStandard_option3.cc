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
/// \file medical/fanoCavity2/src/PhysListEmStandard_option3.cc
/// \brief Implementation of the PhysListEmStandard_option3 class
//
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "PhysListEmStandard_option3.hh"
#include "DetectorConstruction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4UrbanMscModel95.hh"

#include "G4eIonisation.hh"
#include "MyMollerBhabhaModel.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4hIonisation.hh"
#include "G4hMultipleScattering.hh"

#include "G4EmProcessOptions.hh"
#include "G4MscStepLimitType.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard_option3::PhysListEmStandard_option3(const G4String& name,
                               DetectorConstruction* det)
: G4VPhysicsConstructor(name), fDetector(det)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandard_option3::~PhysListEmStandard_option3()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandard_option3::ConstructProcess()
{
  // Add standard EM Processes
  //

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
      // gamma
            
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);
      
    } else if (particleName == "e-") {
      //electron
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->AddEmModel(0, new G4UrbanMscModel95());
            
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetEmModel(new MyMollerBhabhaModel);
                         
      pmanager->AddProcess(msc,                       -1, 1, 1);
      pmanager->AddProcess(eIoni,                     -1, 2, 2);
///      pmanager->AddProcess(new G4eBremsstrahlung,     -1, 3, 3);
            
    } else if (particleName == "e+") {
      //positron
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      msc->AddEmModel(0, new G4UrbanMscModel95());
            
      G4eIonisation* pIoni = new G4eIonisation();
      pIoni->SetEmModel(new MyMollerBhabhaModel);
                               
      pmanager->AddProcess(msc,                       -1, 1, 1);
      pmanager->AddProcess(pIoni,                     -1, 2, 2);
///      pmanager->AddProcess(new G4eBremsstrahlung,     -1, 3, 3);
      pmanager->AddProcess(new G4eplusAnnihilation,    0,-1, 3);
             
    } else if( particleName == "proton" ) {
      //proton  
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
    }
  }

  // Em options
  //
  // Main options and setting parameters are shown here.
  // Several of them have default values.
  //
  G4EmProcessOptions emOptions;
  
  //physics tables
  //
  emOptions.SetMinEnergy(100*eV);        //default    
  emOptions.SetMaxEnergy(10*GeV);        //default  
  emOptions.SetDEDXBinning(8*20);        //default=8*7
  emOptions.SetLambdaBinning(8*20);      //default=8*7
  emOptions.SetSplineFlag(true);         //default
      
  //multiple coulomb scattering
  //
  emOptions.SetMscStepLimitation(fUseDistanceToBoundary);  //default=fUseSafety
      
  //energy loss
  //
  emOptions.SetStepFunction(0.2, 10*um); //default=(0.2, 1*mm)   
           
  //build CSDA range
  //
  emOptions.SetBuildCSDARange(true);          //default=false
  emOptions.SetMaxEnergyForCSDARange(10*GeV); 
  emOptions.SetDEDXBinningForCSDARange(8*20); //default=8*7
          
  //ionization
  //
  emOptions.SetSubCutoff(false);         //default
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

