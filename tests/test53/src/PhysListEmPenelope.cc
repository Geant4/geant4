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
// $Id: PhysListEmPenelope.cc,v 1.2 2009-03-25 13:18:02 pandola Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "PhysListEmPenelope.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4PhotoElectricEffect.hh"
#include "G4PenelopePhotoElectricModel.hh"

#include "G4ComptonScattering.hh"
#include "G4PenelopeComptonModel.hh"

#include "G4GammaConversion.hh"
#include "G4PenelopeGammaConversionModel.hh"

#include "G4RayleighScattering.hh" 
#include "G4PenelopeRayleighModel.hh"

#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4PenelopeIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4PenelopeBremsstrahlungModel.hh"

#include "G4eplusAnnihilation.hh"
#include "G4PenelopeAnnihilationModel.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmPenelope::PhysListEmPenelope(const G4String& name)
  :  G4VPhysicsConstructor(name)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmPenelope::~PhysListEmPenelope()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmPenelope::ConstructProcess()
{
  // Add Penelope or standard EM Processes

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
      // gamma       
      G4PhotoElectricEffect* ePhot = new G4PhotoElectricEffect();
      ePhot->SetModel(new G4PenelopePhotoElectricModel());
      pmanager->AddDiscreteProcess(ePhot);

      G4ComptonScattering* eComp = new G4ComptonScattering();
      eComp->SetModel(new G4PenelopeComptonModel());
      pmanager->AddDiscreteProcess(eComp);

      G4GammaConversion* eGammaConv = new G4GammaConversion();
      eGammaConv->SetModel(new G4PenelopeGammaConversionModel());     
      pmanager->AddDiscreteProcess(eGammaConv);

      G4RayleighScattering* eRayl = new G4RayleighScattering();
      eRayl->SetModel(new G4PenelopeRayleighModel());
      pmanager->AddDiscreteProcess(eRayl);
      
    } else if (particleName == "e-") {
      //electron
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetEmModel(new G4PenelopeIonisationModel());
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
      eBrem->SetEmModel(new G4PenelopeBremsstrahlungModel());
      
      pmanager->AddProcess(new G4MultipleScattering,     -1, 1, 1);
      pmanager->AddProcess(eIoni,     -1, 2, 2);
      pmanager->AddProcess(eBrem, -1,-1, 3);
	    
    } else if (particleName == "e+") {
      //positron
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetEmModel(new G4PenelopeIonisationModel());
      G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
      eBrem->SetEmModel(new G4PenelopeBremsstrahlungModel());
      G4eplusAnnihilation* eAnni = new G4eplusAnnihilation();
      eAnni->SetModel(new G4PenelopeAnnihilationModel());

      pmanager->AddProcess(new G4MultipleScattering,     -1, 1, 1);
      pmanager->AddProcess(eIoni,     -1, 2, 2);
      pmanager->AddProcess(eBrem, -1,-1, 3);      
      pmanager->AddProcess(eAnni,    0,-1, 4);
      
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon  
      pmanager->AddProcess(new G4MultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4MuIonisation,       -1, 2, 2);
      pmanager->AddProcess(new G4MuBremsstrahlung,   -1, 3, 3);
      pmanager->AddProcess(new G4MuPairProduction,   -1, 4, 4);       
     
    } else if( particleName == "alpha" || particleName == "GenericIon" ) { 
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4ionIonisation,      -1, 2,2);
     
    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddProcess(new G4MultipleScattering, -1, 1, 1);
      pmanager->AddProcess(new G4hIonisation,        -1, 2, 2);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

