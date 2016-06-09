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
//
// $Id: PhysListEmG4v52.cc,v 1.3 2006/06/29 16:58:22 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PhysListEmG4v52.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4ComptonScattering52.hh"
#include "G4GammaConversion52.hh"
#include "G4PhotoElectricEffect52.hh"

#include "G4MultipleScattering52.hh"

#include "G4eIonisation52.hh"
#include "G4eBremsstrahlung52.hh"
#include "G4eplusAnnihilation52.hh"

#include "G4MuIonisation52.hh"
#include "G4MuBremsstrahlung52.hh"
#include "G4MuPairProduction52.hh"

#include "G4hIonisation52.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmG4v52::PhysListEmG4v52(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmG4v52::~PhysListEmG4v52()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmG4v52::ConstructProcess()
{
  // Add EM processes realised on base of prototype of model approach design

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
      // gamma         
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect52);
      pmanager->AddDiscreteProcess(new G4ComptonScattering52);
      pmanager->AddDiscreteProcess(new G4GammaConversion52);
      
    } else if (particleName == "e-") {
      //electron
      pmanager->AddProcess(new G4MultipleScattering52, -1, 1,1);
      pmanager->AddProcess(new G4eIonisation52,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung52,    -1,-1,3);
	    
    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4MultipleScattering52, -1, 1,1);
      pmanager->AddProcess(new G4eIonisation52,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung52,    -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation52,   0,-1,4);
      
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon  
      pmanager->AddProcess(new G4MultipleScattering52,-1, 1,1);
      pmanager->AddProcess(new G4MuIonisation52,      -1, 2,2);
      pmanager->AddProcess(new G4MuBremsstrahlung52,  -1,-1,3);
      pmanager->AddProcess(new G4MuPairProduction52,  -1,-1,4);       
     
    } else if( particleName == "GenericIon" ) {
 
      pmanager->AddProcess(new G4MultipleScattering52,-1,1,1);
      pmanager->AddProcess(new G4hIonisation52,       -1,2,2);

    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddProcess(new G4MultipleScattering52,-1,1,1);
      pmanager->AddProcess(new G4hIonisation52,       -1,2,2);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

