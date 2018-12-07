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
/// \file medical/electronScattering/src/PhysListEmStandardSS.cc
/// \brief Implementation of the PhysListEmStandardSS class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysListEmStandardSS.hh"

#include "G4BuilderType.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4CoulombScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandardSS::PhysListEmStandardSS(const G4String& name)
   :  G4VPhysicsConstructor(name)
{
    SetPhysicsType(bElectromagnetic);

    G4EmParameters* param = G4EmParameters::Instance();
    param->SetDefaults();
    param->SetMscThetaLimit(0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandardSS::~PhysListEmStandardSS()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandardSS::ConstructProcess()
{
  // Add standard EM Processes

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
      // gamma
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);
      
    } else if (particleName == "e-") {
      //electron
      pmanager->AddProcess(new G4eIonisation,        -1, 1, 1);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 2, 2);
      pmanager->AddDiscreteProcess(new G4CoulombScattering);            
            
    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4eIonisation,        -1, 1, 1);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 2, 2);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1, 3);
      pmanager->AddDiscreteProcess(new G4CoulombScattering);            
            
    } else if (particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon
      pmanager->AddProcess(new G4MuIonisation,       -1, 1, 1);
      pmanager->AddProcess(new G4MuBremsstrahlung,   -1, 2, 2);
      pmanager->AddProcess(new G4MuPairProduction,   -1, 3, 3);
      pmanager->AddDiscreteProcess(new G4CoulombScattering);              
             
    } else if (particleName == "alpha" || particleName == "He3") {
      G4CoulombScattering* cs = new G4CoulombScattering();
      cs->SetBuildTableFlag(false);
      pmanager->AddProcess(new G4ionIonisation,      -1, 1, 1);      
      pmanager->AddDiscreteProcess(cs);

    } else if (particleName == "GenericIon" ) { 
      G4CoulombScattering* cs = new G4CoulombScattering();
      cs->SetBuildTableFlag(false);
      pmanager->AddProcess(new G4ionIonisation,      -1, 1, 1);      
      pmanager->AddDiscreteProcess(cs);
     
    } else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) && 
               (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddDiscreteProcess(new G4CoulombScattering);            
      pmanager->AddProcess(new G4hIonisation,        -1, 1, 1);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

