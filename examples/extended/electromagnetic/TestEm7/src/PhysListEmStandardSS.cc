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
// $Id: PhysListEmStandardSS.cc,v 1.5 2008/01/14 12:11:39 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01-patch-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "PhysListEmStandardSS.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"
#include "G4ionGasIonisation.hh"
#include "G4IonFluctuations.hh"
#include "G4CoulombScattering.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandardSS::PhysListEmStandardSS(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmStandardSS::~PhysListEmStandardSS()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmStandardSS::ConstructProcess()
{
  // Add standard EM Processes

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
      pmanager->AddProcess(new G4eIonisation,        -1, 1,1);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 2,2);
      pmanager->AddDiscreteProcess(new G4CoulombScattering);      
	    
    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4eIonisation,        -1, 1,1);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1, 2,2);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,3);
      pmanager->AddDiscreteProcess(new G4CoulombScattering);      
            
    } else if (particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon  
      pmanager->AddProcess(new G4MuIonisation,       -1, 1,1);
      pmanager->AddProcess(new G4MuBremsstrahlung,   -1, 2,2);
      pmanager->AddProcess(new G4MuPairProduction,   -1, 3,3);
      pmanager->AddDiscreteProcess(new G4CoulombScattering);      
             
    } else if (particleName == "alpha" || particleName == "He3") {
      G4ionIonisation* ion = new G4ionIonisation();
      ion->SetStepFunction(0.1, um);
      ion->ActivateNuclearStopping(false);
      pmanager->AddProcess(ion,  -1, 1,1);
      pmanager->AddDiscreteProcess(new G4CoulombScattering);      

    } else if (particleName == "GenericIon" ) { 
      G4ionGasIonisation* ion = new G4ionGasIonisation();
      ion->ActivateNuclearStopping(false);
      ion->SetStepFunction(0.1, um);
      pmanager->AddProcess(ion,  -1, 1,1);
      G4CoulombScattering* cs = new G4CoulombScattering();
      cs->SetBuildTableFlag(false);
      pmanager->AddDiscreteProcess(cs);
      
    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      G4hIonisation* hion = new G4hIonisation();
      hion->SetStepFunction(0.1, 10.*um);
      hion->ActivateNuclearStopping(false);
      pmanager->AddProcess(hion,  -1,1,1);
      pmanager->AddDiscreteProcess(new G4CoulombScattering);      
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

