//---------------------------------------------------------------------------
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
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

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCaptureAtRest.hh"

#include "G4hIonisation.hh"
#include "hTestStepCut.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestStanPhysicsList::hTestStanPhysicsList():
  verbose(0),
  maxChargedStep(DBL_MAX)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestStanPhysicsList::ConstructProcess()
{
  hTestStepCut* theStepCut = new hTestStepCut();
  theStepCut->SetMaxStep(maxChargedStep);          

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {

      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      pmanager->AddDiscreteProcess(new G4ComptonScattering);
      pmanager->AddDiscreteProcess(new G4GammaConversion);    
      
    } else if (particleName == "e-") {
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);

      pmanager->AddProcess(new G4eIonisation,   -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);       
      pmanager->AddProcess(theStepCut,  -1,-1,4);

    } else if (particleName == "e+") {

      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);
      pmanager->AddProcess(theStepCut,   -1,-1,5);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {

      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4MuIonisation,      -1, 2,2);
      pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1,3);
      pmanager->AddProcess(new G4MuPairProduction,  -1,-1,4);       	       
      pmanager->AddProcess(new G4MuonMinusCaptureAtRest,0,-1,-1);
      pmanager->AddProcess(theStepCut,   -1,-1,5);    

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

      pmanager->AddProcess(hIon,-1,2,2);
      pmanager->AddProcess(theStepCut,  -1,-1,3);
   
    } else if (   particleName == "alpha"  
               || particleName == "deuteron"  
               || particleName == "triton"  
               || particleName == "IonC12"  
               || particleName == "IonAr40"  
               || particleName == "IonFe56"  
              )
    {
      pmanager->AddProcess(new G4MultipleScattering,-1,1,1);

      if(0 < verbose) {
        G4cout << "Ionic processes for " << particleName << G4endl; 
      }

      G4hIonisation* iIon = new G4hIonisation() ;

      pmanager->AddProcess(iIon,-1,2,2);
      pmanager->AddProcess(theStepCut,  -1,-1,3);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


