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
// $Id: PhysListEmPAI.cc,v 1.1 2004/05/27 18:04:31 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-02 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#include "PhysListEmPAI.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4PAIModel.hh"
#include "G4PAIPhotonModel.hh"
//#include "G4PAIwithPhotons.hh"

#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmPAI::PhysListEmPAI(const G4String& name)
   :  G4VPhysicsConstructor(name)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysListEmPAI::~PhysListEmPAI()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysListEmPAI::ConstructProcess()
{
  // Add standard EM Processes

  const G4RegionStore* theRegionStore = G4RegionStore::GetInstance();
  G4Region* gas = theRegionStore->GetRegion("VertexDetector");

  theParticleIterator->reset();

  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager     = particle->GetProcessManager();
    G4String particleName          = particle->GetParticleName();
     
    if (particleName == "gamma") 
    {     
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);
      
    } 
    else if (particleName == "e-") 
    { 
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);

      G4eIonisation* eion = new G4eIonisation();


      G4PAIModel*     pai = new G4PAIModel(particle,"PAIModel");
      //G4PAIPhotonModel*     pai = new G4PAIPhotonModel(particle,"PAIModel");
      // G4PAIwithPhotons*     pai   = new G4PAIwithPhotons(particle,"PAIwPhotons");
      eion->AddEmModel(0,pai,pai,gas);

      pmanager->AddProcess(eion,-1, 2,2);

      pmanager->AddProcess(new G4eBremsstrahlung,-1,1,3);
	    
    } 
    else if (particleName == "e+") 
    {
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);
      
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) 
    {
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);

      G4MuIonisation* muion = new G4MuIonisation();

      G4PAIModel*     pai   = new G4PAIModel(particle,"PAIModel");
      // G4PAIwithPhotons*     pai   = new G4PAIwithPhotons(particle,"PAIwPhotons");
      muion->AddEmModel(0,pai,pai,gas);

      pmanager->AddProcess(muion,      -1, 2,-2);
      pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1,3);
      pmanager->AddProcess(new G4MuPairProduction,  -1,-1,4);       

    } 
    else if (particleName == "GenericIon") 
    {
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4ionIonisation,      -1, 2,2);
     
    } 
    else if ( particleName == "proton" ) 
    {

      pmanager->AddProcess(new G4MultipleScattering,-1,1,1);

      G4hIonisation* pion =     new G4hIonisation();
      G4PAIModel*     pai = new G4PAIModel(particle,"PAIModel");
      //G4PAIPhotonModel*     pai = new G4PAIPhotonModel(particle,"PAIModel");
      // G4PAIwithPhotons*     pai   = new G4PAIwithPhotons(particle,"PAIwPhotons");
      pion->AddEmModel(0,pai,pai,gas);

      pmanager->AddProcess(pion,       -1,2,2);
    }
    else if ( ( !particle->IsShortLived() )       &&
	      ( particle->GetPDGCharge() != 0.0 ) && 
	      ( particle->GetParticleName() != "chargedgeantino") ) 
    {

      pmanager->AddProcess(new G4MultipleScattering,-1,1,1);

      G4hIonisation* hion =     new G4hIonisation();

      // G4PAIModel*     pai = new G4PAIModel(particle,"PAIModel");
      //  G4PAIwithPhotons*     pai   = new G4PAIwithPhotons(particle,"PAIwPhotons");
      // hion->AddEmModel(0,pai,pai,gas);

      pmanager->AddProcess(hion,       -1,2,2);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

