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
//---------------------------------------------------------------------------
//
//---------------------------------------------------------------------------
//
// ClassName:   hTestStanPhysicsList
//  
// Description: Standard EM physics list
//
// Authors:   08.04.01  V.Ivanchenko 
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestStanPhysicsList.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScattering.hh"
#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCaptureAtRest.hh"

#include "G4hIonisation.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestStanPhysicsList::hTestStanPhysicsList()
{
  verbose = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestStanPhysicsList::ConstructProcess()
{
  G4cout << "Standard EM PhysicsList will be constructed" << G4endl;

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    G4String particleType = particle->GetParticleType();
    G4double charge = particle->GetPDGCharge();
     
    if (particleName == "gamma") {

      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);    
      
    } else if (particleName == "e-") {
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      G4eIonisation* eion = new G4eIonisation();
      eion->SetSubSec(false);
      pmanager->AddProcess(eion,   -1, 2,2);
      G4eBremsstrahlung* ebr = new G4eBremsstrahlung();
      ebr->SetSubSec(false);
      pmanager->AddProcess(ebr,    -1,-1,3);       

    } else if (particleName == "e+") {

      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      G4eIonisation* pion = new G4eIonisation();
      pion->SetSubSec(false);
      pmanager->AddProcess(pion,        -1, 2,2);
      G4eBremsstrahlung* pbr = new G4eBremsstrahlung();
      pbr->SetSubSec(false);
      pmanager->AddProcess(pbr,    -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {

      
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4MuIonisation,      -1, 2,2);
      pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1,3);
      pmanager->AddProcess(new G4MuPairProduction,  -1,-1,4);       	       
      pmanager->AddProcess(new G4MuonMinusCaptureAtRest,0,-1,-1);
      
    } else if (
                particleName == "proton"  
               || particleName == "anti_proton"  
               || particleName == "pi+"  
               || particleName == "pi-"  
               || particleName == "kaon+"  
               || particleName == "kaon-"  
              )
    {
      pmanager->AddProcess(new G4MultipleScattering,-1,1,1);
      if(0 < verbose) {
        G4cout << "Hadronic processes for " << particleName << G4endl; 
      }

      G4hIonisation* hIon = new G4hIonisation() ;
      hIon->SetVerboseLevel(0);
      pmanager->AddProcess(hIon,-1,2,2);
   
    } else if (   particleName == "alpha"  
               || particleName == "deuteron"  
               || particleName == "triton"  
               || particleName == "IonC12"  
               || particleName == "IonAr40"  
               || particleName == "IonFe56"  
               || particleName == "GenericIon"
               || particleType == "nucleus") 
              
    {
      pmanager->AddProcess(new G4MultipleScattering,-1,1,1);

      if(0 < verbose) {
        G4cout << "Ionic processes for " << particleName << G4endl; 
      }

      G4hIonisation* iIon = new G4hIonisation();
      iIon->SetVerboseLevel(1);
      pmanager->AddProcess(iIon,-1,1,1);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


