//---------------------------------------------------------------------------
//
// ClassName:   hTestLowEPhysicsList
//  
// Description: LowEnergy EM processes list
//
// Authors:     V.Ivanchenko 29/03/01
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

void hTestLowEPhysicsList::hTestLowEPhysicsList():
  verbose(0),
  maxChargedStep(DBL_MAX),
  nuclStop(true),
  barkas(true),
  table("")
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestLowEPhysicsList::ConstructProcess()
{
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

      hTestStepCut* theeminusStepCut = new hTestStepCut();
      theeminusStepCut->SetMaxStep(maxChargedStep);  
      pmanager->AddProcess(theeminusStepCut,         -1,-1,4);

    } else if (particleName == "e+") {
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      pmanager->AddProcess(new G4LowEnergyIonisation,  -1, 2,2);
      pmanager->AddProcess(new G4LowEnergyBremsstrahlung, -1,-1,3);   
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);
                  
      hTestStepCut* theeplusStepCut = new hTestStepCut();
      theeplusStepCut->SetMaxStep(maxChargedStep) ;          
      pmanager->AddProcess(theeplusStepCut,          -1,-1,5);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
      pmanager->AddProcess(new G4MultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4MuIonisation,      -1, 2,2);
      pmanager->AddProcess(new G4MuBremsstrahlung,  -1,-1,3);
      pmanager->AddProcess(new G4MuPairProduction,  -1,-1,4);       	       
      pmanager->AddProcess(new G4MuonMinusCaptureAtRest,0,-1,-1);
	       
      hTestStepCut* muStepCut = new hTestStepCut();
      muStepCut->SetMaxStep(maxChargedStep) ;          
      pmanager->AddProcess(muStepCut,          -1,-1,5);

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

      G4hLowEnergyIonisation* hIon = new G4hLowEnergyIonisation() ;

      if(nuclStop) hIon->SetNuclearStoppingOn();
      else         hIon->SetNuclearStoppingOff();

      if(barkas)   hIon->SetBarkasOn();
      else         hIon->SetBarkasOff();

      if(table == G4String("Ziegler1977He") ||
         table == G4String("Ziegler1977H") ||
         table == G4String("ICRU_R49p") ||
         table == G4String("ICRU_R49He") ||
         table == G4String("ICRU_R49PowersHe") )
         hIon->SetStoppingPowerTableName(table);

      pmanager->AddProcess(hIon,-1,2,2);
      
      hTestStepCut* thehadronStepCut = new hTestStepCut();
      thehadronStepCut->SetMaxStep(maxChargedStep);          		       
      pmanager->AddProcess( thehadronStepCut,       -1,-1,3);
   
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
         iIon->SetStoppingPowerTableName(table);

      pmanager->AddProcess(iIon,-1,2,2);
      
      hTestStepCut* theIonStepCut = new hTestStepCut();
      theIonStepCut->SetMaxStep(maxChargedStep);     
      pmanager->AddProcess( theIonStepCut,       -1,-1,3);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


