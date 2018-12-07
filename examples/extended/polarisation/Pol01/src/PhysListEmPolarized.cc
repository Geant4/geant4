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
/// \file polarisation/Pol01/src/PhysListEmPolarized.cc
/// \brief Implementation of the PhysListEmPolarized class
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysListEmPolarized.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4eMultipleScattering.hh"

#include "G4PolarizedCompton.hh"
#include "G4PolarizedGammaConversion.hh"
#include "G4ePolarizedIonisation.hh"
#include "G4ePolarizedBremsstrahlung.hh"
#include "G4eplusPolarizedAnnihilation.hh"
#include "G4PolarizedPhotoElectricEffect.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmPolarized::PhysListEmPolarized(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmPolarized::~PhysListEmPolarized()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmPolarized::ConstructProcess()
{
  // Add standard EM Processes
  
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();

    if (particleName == "gamma") {
      pmanager->AddDiscreteProcess(new G4PolarizedPhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4PolarizedCompton);
      pmanager->AddDiscreteProcess(new G4PolarizedGammaConversion);      

    } else if (particleName == "e-") {
      pmanager->AddProcess(new G4eMultipleScattering,   -1,1,1);
      pmanager->AddProcess(new G4ePolarizedIonisation,  -1,2,2);
      pmanager->AddProcess(new G4ePolarizedBremsstrahlung,      -1,3,3);

    } else if (particleName == "e+") {
      pmanager->AddProcess(new G4eMultipleScattering,  -1, 1,1);
      pmanager->AddProcess(new G4ePolarizedIonisation, -1, 2,2);
      pmanager->AddProcess(new G4ePolarizedBremsstrahlung,    -1, 3,3);
      pmanager->AddProcess(new G4eplusPolarizedAnnihilation,   0,-1,4);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
