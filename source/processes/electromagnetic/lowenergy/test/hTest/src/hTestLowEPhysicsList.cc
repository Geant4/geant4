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
// ClassName:   hTestLowEPhysicsList
//  
// Description: LowEnergy EM processes list
//
// Authors:    08.04.01 V.Ivanchenko 
//
// Modified:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "hTestLowEPhysicsList.hh"

#include "G4LowEnergyRayleigh.hh"
#include "G4LowEnergyCompton.hh"   
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"

#include "G4MultipleScattering.hh"

#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCaptureAtRest.hh"

#include "G4hLowEnergyIonisation.hh"
#include "hTestStepCut.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestLowEPhysicsList::hTestLowEPhysicsList()
{
  verbose = 0;
  maxChargedStep = DBL_MAX;
  nuclStop = true;
  barkas = true;
  table = G4String("");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestLowEPhysicsList::ConstructProcess()
{
  G4cout << "LowEnergy Electromagnetic PhysicsList is initilized" << G4endl;
  hTestStepCut* theStepCut = new hTestStepCut();
  theStepCut->SetMaxStep(maxChargedStep);
  G4ParticleDefinition* proton = G4Proton::ProtonDefinition();     
  G4ParticleDefinition* antiproton = G4AntiProton::AntiProtonDefinition();     

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
      pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh);
      pmanager->AddDiscreteProcess(new G4LowEnergyPhotoElectric);
      pmanager->AddDiscreteProcess(new G4LowEnergyCompton);
      pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion);    
      
    } else if (particleName == "e-") {
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4LowEnergyIonisation,  -1, 2,2);
      pmanager->AddProcess(new G4LowEnergyBremsstrahlung, -1,-1,3);   
      pmanager->AddProcess(theStepCut, -1,-1,4);

    } else if (particleName == "e+") {
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4LowEnergyIonisation,  -1, 2,2);
      pmanager->AddProcess(new G4LowEnergyBremsstrahlung, -1,-1,3);   
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);
      pmanager->AddProcess(theStepCut, -1,-1,5);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4MuIonisation,      -1, 2,2);
      pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1,3);
      pmanager->AddProcess(new G4MuPairProduction,  -1,-1,4);       	       
      pmanager->AddProcess(new G4MuonMinusCaptureAtRest,0,-1,-1);
      pmanager->AddProcess(theStepCut, -1,-1,5);

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
        G4cout << "Hadronic EM processes for " << particleName << G4endl; 
      }

      G4hLowEnergyIonisation* hIon = new G4hLowEnergyIonisation() ;

      if(nuclStop) hIon->SetNuclearStoppingOn();
      else         hIon->SetNuclearStoppingOff();

      if(barkas)   hIon->SetBarkasOn();
      else         hIon->SetBarkasOff();

      if(table == G4String("Ziegler1977He") ||
         table == G4String("Ziegler1977H") ||
         table == G4String("ICRU_R49p") ||
         table == G4String("ICRU_R49He") ||
         table == G4String("ICRU_R49PowersHe") ) {

	if(particle->GetPDGCharge() > 0.0) 
          hIon->SetElectronicStoppingPowerModel(proton,table);
        else
          hIon->SetElectronicStoppingPowerModel(antiproton,table);
      }

      pmanager->AddProcess(hIon,-1,2,2);		       
      pmanager->AddProcess( theStepCut,  -1,-1,3);
   
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
        G4cout << "Ionic EM processes for " << particleName << G4endl; 
      }

      G4hLowEnergyIonisation* iIon = new G4hLowEnergyIonisation() ;

      if(nuclStop) iIon->SetNuclearStoppingOn();
      else         iIon->SetNuclearStoppingOff();

      if(barkas)   iIon->SetBarkasOn();
      else         iIon->SetBarkasOff();

      if(table == G4String("Ziegler1977He") ||
         table == G4String("Ziegler1977H") ||
         table == G4String("ICRU_R49p") ||
         table == G4String("ICRU_R49He") ||
         table == G4String("ICRU_R49PowersHe") )
         iIon->SetElectronicStoppingPowerModel(proton,table);

      pmanager->AddProcess(iIon,-1,2,2);
      pmanager->AddProcess(theStepCut, -1,-1,3);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


