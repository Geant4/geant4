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
/// \file electromagnetic/TestEm5/src/PhysListEmStandardSSM.cc
/// \brief Implementation of the PhysListEmStandardSSM class
//
// $Id: PhysListEmStandardSSM.cc 100286 2016-10-17 08:43:45Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysListEmStandardSSM.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4CoulombScattering.hh"
#include "G4eCoulombScatteringModel.hh"
#include "G4hCoulombScatteringModel.hh"
#include "G4eSingleCoulombScatteringModel.hh"
#include "G4IonCoulombScatteringModel.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandardSSM::PhysListEmStandardSSM(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandardSSM::~PhysListEmStandardSSM()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandardSSM::ConstructProcess()
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

      G4CoulombScattering* cs = new G4CoulombScattering();
      G4eSingleCoulombScatteringModel* model = 
        new G4eSingleCoulombScatteringModel();
      //model->SetLowEnergyThreshold(10*eV);
      model->SetPolarAngleLimit(0.0);
      cs->AddEmModel(0, model);
      pmanager->AddDiscreteProcess(cs);            
            
    } else if (particleName == "e+") {
      //positron

      pmanager->AddProcess(new G4eIonisation,        -1, 1, 1);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 2, 2);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1, 3);

      G4CoulombScattering* cs = new G4CoulombScattering();
      G4eSingleCoulombScatteringModel* model = 
        new G4eSingleCoulombScatteringModel();
      model->SetPolarAngleLimit(0.0);
      cs->AddEmModel(0, model);
      pmanager->AddDiscreteProcess(cs);            
            
    } else if (particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon
      pmanager->AddProcess(new G4MuIonisation,       -1, 1, 1);
      pmanager->AddProcess(new G4MuBremsstrahlung,   -1, 2, 2);
      pmanager->AddProcess(new G4MuPairProduction,   -1, 3, 3);
      G4CoulombScattering* cs = new G4CoulombScattering();
      G4hCoulombScatteringModel* model = new G4hCoulombScatteringModel();
      model->SetPolarAngleLimit(0.0);
      cs->AddEmModel(0, model);
      pmanager->AddDiscreteProcess(cs);            

             
    } else if (particleName == "alpha" || particleName == "He3") {
      pmanager->AddProcess(new G4ionIonisation,      -1, 1, 1);
      G4CoulombScattering* cs = new G4CoulombScattering();
      cs->AddEmModel(0, new G4IonCoulombScatteringModel());
      cs->SetBuildTableFlag(false);
      pmanager->AddDiscreteProcess(cs);            

    } else if (particleName == "GenericIon" ) { 
      pmanager->AddProcess(new G4ionIonisation,      -1, 1, 1);      
      G4CoulombScattering* cs = new G4CoulombScattering();
      cs->AddEmModel(0, new G4IonCoulombScatteringModel());
      cs->SetBuildTableFlag(false);
      pmanager->AddDiscreteProcess(cs);
     
    } else if ((!particle->IsShortLived()) &&
               (particle->GetPDGCharge() != 0.0) && 
               (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddProcess(new G4hIonisation,        -1, 1, 1);
      pmanager->AddDiscreteProcess(new G4CoulombScattering);            
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

