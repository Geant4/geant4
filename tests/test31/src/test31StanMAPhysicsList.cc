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
// ClassName:   test31StanMAPhysicsList
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

#include "test31StanMAPhysicsList.hh"

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4MultipleScatteringSTD.hh"

#include "G4eIonisationSTD.hh"
#include "G4eBremsstrahlungSTD.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisationSTD.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCaptureAtRest.hh"

#include "G4hIonisationSTD.hh"
#include "G4ionIonisation.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

test31StanMAPhysicsList::test31StanMAPhysicsList()
{
  verbose = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void test31StanMAPhysicsList::ConstructProcess()
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
      //      pmanager->AddProcess(new G4MultipleScatteringSTD, -1, 1,1);
      G4eIonisationSTD* eion = new G4eIonisationSTD();
      eion->SetLossFluctuations(false);
      eion->SetVerboseLevel(1);
      //      eion->SetParticle(particle);
      //      eion->SetSubCutoff(false);
      pmanager->AddProcess(eion,   -1, 2,2);
      G4eBremsstrahlungSTD* ebr = new G4eBremsstrahlungSTD();
      //      ebr->SetParticle(particle);
      //      ebr->SetSubCutoff(false);
      pmanager->AddProcess(ebr,    -1,-1,3);       

    } else if (particleName == "e+") {

      pmanager->AddProcess(new G4MultipleScatteringSTD, -1, 1,1);
      G4eIonisationSTD* pion = new G4eIonisationSTD();
      //      pion->SetParticle(particle);
      //      pion->SetSubCutoff(false);
      pmanager->AddProcess(pion,        -1, 2,2);
      G4eBremsstrahlungSTD* pbr = new G4eBremsstrahlungSTD();
      //      pbr->SetParticle(particle);
      //      pbr->SetSubCutoff(false);
      pmanager->AddProcess(pbr,    -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {

      pmanager->AddProcess(new G4MultipleScatteringSTD,-1, 1,1);
      pmanager->AddProcess(new G4MuIonisationSTD,      -1, 2,2);
      //      pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1,3);
      //      pmanager->AddProcess(new G4MuPairProduction,  -1,-1,4);       	       
      //      pmanager->AddProcess(new G4MuonMinusCaptureAtRest,0,-1,-1);

    } else if (
                particleName == "proton"  
               || particleName == "anti_proton"  
              )
    {
      pmanager->AddProcess(new G4MultipleScatteringSTD,-1,1,1);
      if(0 < verbose) {
        G4cout << "Hadronic processes for " << particleName << G4endl; 
      }

      G4hIonisationSTD* hion = new G4hIonisationSTD();
      //      hion->SetParticle(particle);
      //      hion->SetSubCutoff(false);
      pmanager->AddProcess(hion,-1,2,2);

   } else if (
                particleName == "pi+"  
               || particleName == "kaon+"   
              )
    {
      pmanager->AddProcess(new G4MultipleScatteringSTD,-1,1,1);
      if(0 < verbose) {
        G4cout << "Hadronic processes for " << particleName << G4endl; 
      }

      G4hIonisationSTD* hion = new G4hIonisationSTD();
      //      hion->SetParticle(particle);
      //      hion->SetSubCutoff(false);
      //      hion->SetBaseParticle(G4Proton::Proton());
      pmanager->AddProcess(hion,-1,2,2);   

   } else if (
                particleName == "pi-"      
               || particleName == "kaon-"  
              )
    {
      pmanager->AddProcess(new G4MultipleScatteringSTD,-1,1,1);
      if(0 < verbose) {
        G4cout << "Hadronic processes for " << particleName << G4endl; 
      }

      G4hIonisationSTD* hion = new G4hIonisationSTD();
      //      hion->SetParticle(particle);
      //      hion->SetSubCutoff(false);
      //      hion->SetBaseParticle(G4AntiProton::AntiProton());
      pmanager->AddProcess(hion,-1,2,2);   


    } else if (   particleName == "alpha"  
               || particleName == "deuteron"  
               || particleName == "triton"  
               || particleName == "IonC12"  
               || particleName == "IonAr40"  
               || particleName == "IonFe56")
    {
      pmanager->AddProcess(new G4MultipleScatteringSTD,-1,1,1);

      if(0 < verbose) {
        G4cout << "Ionic processes for " << particleName << G4endl; 
      }

      G4hIonisationSTD* iIon = new G4hIonisationSTD();
      //      iIon->SetParticle(particle);
      //      iIon->SetSubCutoff(false);
      //      iIon->SetBaseParticle(G4Proton::Proton());
      pmanager->AddProcess(iIon,-1,2,2);

    } else if (particleName == "GenericIon"  
               || (particleType == "nucleus" && charge != 0.0) )
    {
      pmanager->AddProcess(new G4MultipleScatteringSTD,-1,1,1);

      if(0 < verbose) {
        G4cout << "Ionic processes for " << particleName << G4endl; 
      }

      G4ionIonisation* iIon = new G4ionIonisation();
      //      iIon->SetParticle(particle);
      //      iIon->SetSubCutoff(false);
      //      iIon->SetBaseParticle(G4Proton::Proton());
      pmanager->AddProcess(iIon,-1,2,2);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


