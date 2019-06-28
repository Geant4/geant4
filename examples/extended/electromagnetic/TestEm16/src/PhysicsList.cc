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
/// \file electromagnetic/TestEm16/src/PhysicsList.cc
/// \brief Implementation of the PhysicsList class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

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

#include "G4SynchrotronRadiation.hh"
#include "G4SynchrotronRadiationInMat.hh"

#include "G4StepLimiter.hh"
#include "G4DecayPhysics.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList()
: G4VUserPhysicsList(),
  fMess(0)
{
  SetDefaultCutValue(1.*km);
  
  fSRType = true; 
  fMess = new PhysicsListMessenger(this);
  fDecayPhysics = new G4DecayPhysics();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList()
{ 
  delete fMess;
  delete fDecayPhysics;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
  fDecayPhysics->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  fDecayPhysics->ConstructProcess();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructEM()
{
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  G4SynchrotronRadiation* fSync = new G4SynchrotronRadiation();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      // gamma
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);
      pmanager->AddDiscreteProcess(new G4RayleighScattering);
            
    } else if (particleName == "e-") {
      //electron
      pmanager->AddProcess(new G4eMultipleScattering,       -1, 1, -1);
      pmanager->AddProcess(new G4eIonisation,               -1, 2, 1);
      pmanager->AddProcess(new G4eBremsstrahlung,           -1, 3, 2);
      if (fSRType) {
        pmanager->AddProcess(fSync,    -1,-1, 3);
      } else {
        pmanager->AddProcess(new G4SynchrotronRadiationInMat,-1,-1, 3); 
      }
      pmanager->AddProcess(new G4StepLimiter,               -1,-1, 4);
     
    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4eMultipleScattering,       -1, 1, -1);
      pmanager->AddProcess(new G4eIonisation,               -1, 2, 1);
      pmanager->AddProcess(new G4eBremsstrahlung,           -1, 3, 2);
      pmanager->AddProcess(new G4eplusAnnihilation,          0,-1, 3);
      if (fSRType) {
        pmanager->AddProcess(fSync,      -1,-1, 4);
      } else {
        pmanager->AddProcess(new G4SynchrotronRadiationInMat, -1,-1, 4);
      }
      pmanager->AddProcess(new G4StepLimiter,               -1,-1, 5);
      
    } else if( particleName == "mu+" ||
               particleName == "mu-"    ) {
      //muon
      pmanager->AddProcess(new G4MuMultipleScattering, -1, 1, -1);
      pmanager->AddProcess(new G4MuIonisation,         -1, 2, 1);
      pmanager->AddProcess(new G4MuBremsstrahlung,     -1, 3, 2);
      pmanager->AddProcess(new G4MuPairProduction,     -1, 4, 3);
      pmanager->AddProcess(fSync, -1,-1, 4); 
      pmanager->AddProcess(new G4StepLimiter,               -1,-1, 5);

    } else if( particleName == "proton") {
      //proton
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, -1);
      pmanager->AddProcess(new G4hIonisation,         -1, 2, 1);
      pmanager->AddProcess(new G4hBremsstrahlung,     -1, 3, 2);
      pmanager->AddProcess(new G4hPairProduction,     -1, 4, 3);
      pmanager->AddProcess(fSync, -1,-1, 4); 
      pmanager->AddProcess(new G4StepLimiter,          -1,-1, 5);
    }
    else if (particle->GetPDGCharge() != 0.0 && !particle->IsShortLived()) 
    {
      pmanager->AddProcess(fSync,-1,-1, 1);  
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
