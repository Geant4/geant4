//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: Em2PhysListEmModel.cc,v 1.2 2002-12-16 11:29:23 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em2PhysListEmModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScatteringSTD.hh"

#include "G4eIonisationSTD.hh"
#include "G4eBremsstrahlungSTD.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisationSTD.hh"
#include "G4MuBremsstrahlungSTD.hh"
#include "G4MuPairProductionSTD.hh"

#include "G4hIonisationSTD.hh"
#include "G4ionIonisation.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em2PhysListEmModel::Em2PhysListEmModel(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em2PhysListEmModel::~Em2PhysListEmModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em2PhysListEmModel::ConstructProcess()
{
  // Add EM processes realised on base of prototype of model approach design

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
      pmanager->AddProcess(new G4MultipleScatteringSTD, -1, 1,1);
      pmanager->AddProcess(new G4eIonisationSTD,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlungSTD,    -1,-1,3);
	    
    } else if (particleName == "e+") {
      //positron
      pmanager->AddProcess(new G4MultipleScatteringSTD, -1, 1,1);
      pmanager->AddProcess(new G4eIonisationSTD,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlungSTD,    -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,      0,-1,4);
      
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      //muon  
      pmanager->AddProcess(new G4MultipleScatteringSTD,-1, 1,1);
      pmanager->AddProcess(new G4MuIonisationSTD,      -1, 2,2);
      pmanager->AddProcess(new G4MuBremsstrahlungSTD,  -1,-1,3);
      pmanager->AddProcess(new G4MuPairProductionSTD,  -1,-1,4);       
     
    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
      //all others charged particles except geantino
      pmanager->AddProcess(new G4MultipleScatteringSTD,-1,1,1);
      if(particle->GetPDGMass() < 1.0*GeV ) {
        pmanager->AddProcess(new G4hIonisationSTD,     -1,2,2);
      } else {
        pmanager->AddProcess(new G4ionIonisation,      -1,2,2);
      }           
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

