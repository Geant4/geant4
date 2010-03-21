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
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

//
// Always place it in front!
//
#include "G4Timer.hh"

#include "Test17PhysicsList.hh"
#include "Test17DetectorConstruction.hh"
#include "Test17PhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProductionCutsTable.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
//#include "G4Material.hh"
//#include "G4EnergyLossTables.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17PhysicsList::Test17PhysicsList(Test17DetectorConstruction* p)
:G4VUserPhysicsList()
{
  pDet = p;

  defaultCutValue = 0.1*mm;
  cutForGamma     = defaultCutValue;
  cutForElectron  = defaultCutValue;

  //  MaxChargedStep = DBL_MAX; 
  MaxChargedStep = 0.1*mm; 

  SetVerboseLevel(2);
  physicsListMessenger = new Test17PhysicsListMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test17PhysicsList::~Test17PhysicsList()
{
  delete physicsListMessenger; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  G4BosonConstructor  pBosonConstructor;
  pBosonConstructor.ConstructParticle();

  G4LeptonConstructor pLeptonConstructor;
  pLeptonConstructor.ConstructParticle();

  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle(); 

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4LowEnergyBremsstrahlung.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
#include "G4hLowEnergyIonisation.hh"

#include "Test17StepCut.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  theStepCut = new Test17StepCut();
  theStepCut->SetMaxStep(MaxChargedStep);
          		       
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
      pmanager->AddProcess(new G4eMultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);       

      pmanager->AddProcess(theStepCut,         -1,-1,4);

    } else if (particleName == "e+") {
      //positron      
      pmanager->AddProcess(new G4eMultipleScattering,-1, 1,1);
      pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);

      pmanager->AddProcess(theStepCut,          -1,-1,5);

    } else if( particleName == "mu+" ||
               particleName == "mu-"    ) {
     //muon
     pmanager->AddProcess(new G4hMultipleScattering,-1, 1,1);
     pmanager->AddProcess(new G4MuIonisation,       -1, 2,2);
     pmanager->AddProcess(new G4MuBremsstrahlung,   -1,-1,3);
     pmanager->AddProcess(new G4MuPairProduction,   -1,-1,4);

    } else if (
                particleName == "proton"  
               || particleName == "anti_proton"  
               || particleName == "pi+"  
               || particleName == "pi-"  
               || particleName == "kaon+"  
               || particleName == "kaon-"  
              )
    {
      pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);

      G4cout << "Hadronic processes for " << particleName << G4endl; 

      // Standard ionisation
      //       G4hIonisation* hIon = new G4hIonisation() ;

      // Standard ionisation with low energy extention
      G4hLowEnergyIonisation* hIon = new G4hLowEnergyIonisation() ;
      //      hIon->SetNuclearStoppingOff() ;
      //         hIon->SetBarkasOff() 
      // hIon->SetNuclearStoppingOff() ;
      // hIon->SetAntiProtonStoppingOff() ;
      //      hIon->SetBarkasOff() ;
      //      hIon->SetNuclearStoppingOn() ;
      hIon->SetFluorescence(true) ;
      hIon->SetCutForSecondaryPhotons(100.*eV) ;

      //hIon->SetStoppingPowerTableName("Ziegler1977He") ;
     //   hIon->SetElectronicStoppingPowerModel(particle,"Ziegler1977p") ;
      //       hIon->SetNuclearStoppingPowerModel("ICRU_R49Mollere") ;
      //   hIon->SetElectronicStoppingPowerModel(particle,"ICRU_R49p") ;

      pmanager->AddProcess(hIon,-1,2,2);
      hionVector.push_back(hIon);
      
      pmanager->AddProcess( theStepCut,       -1,-1,3);

    } else if (   particleName == "alpha"  
		  || particleName == "deuteron"  
		  || particleName == "GenericIon"  
		  )
    {
      pmanager->AddProcess(new G4hMultipleScattering,-1,1,1);

      G4cout << "Ionic processes for " << particleName << G4endl; 

      // Standard ionisation
      //      G4hIonisation* iIon = new G4hIonisation() ;

      // Standard ionisation with low energy extantion
       G4hLowEnergyIonisation* iIon = new G4hLowEnergyIonisation() ;
       iIon->SetVerboseLevel(1);
	   //   iIon->SetNuclearStoppingOff() ;
	//  iIon->SetNuclearStoppingOn() ;

      //iIon->SetStoppingPowerTableName("Ziegler1977He") ;
      //iIon->SetStoppingPowerTableName("Ziegler1977H") ;
	// iIon->SetStoppingPowerTableName("ICRU_R49p") ;
      //iIon->SetStoppingPowerTableName("ICRU_R49He") ;
      //iIon->SetStoppingPowerTableName("ICRU_R49PowersHe") ;
       iIon->SetCutForSecondaryPhotons(100.*eV) ;

       pmanager->AddProcess(iIon,-1,2,2);
       hionVector.push_back(iIon);
      
       pmanager->AddProcess( theStepCut,       -1,-1,3);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4Decay.hh"

void Test17PhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17PhysicsList::SetCuts()
{

  //special for low energy physics
  //
  G4double lowlimit=1*eV;
  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(lowlimit,100*GeV);

  if (verboseLevel >0){
    G4cout << "Test17PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(MaxChargedStep,"Length") << G4endl;
  }
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma

  G4cout << "Set cuts for all particles! " << G4endl;

  SetCutValue(cutForGamma,"gamma");
  SetCutValue(cutForElectron,"e-");
  SetCutValue(cutForElectron,"e+");
  SetCutValue(cutForElectron,"proton");

  if (verboseLevel>0) DumpCutValuesTable();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17PhysicsList::SetGammaCut(G4double val)
{
  cutForGamma = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17PhysicsList::SetElectronCut(G4double val)
{
  cutForElectron = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17PhysicsList::SetMaxStep(G4double step)
{
  MaxChargedStep = step ;
  theStepCut->SetMaxStep(MaxChargedStep);
  G4cout << "Test17PhysicsList:MaxChargedStep=" << MaxChargedStep << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17PhysicsList::SetCutForSecondaryPhotons(G4double cut)
{
  size_t n = hionVector.size();
  G4cout << " Cut for secondary gamma = " << cut/keV << " keV" << G4endl;

  for (size_t i=0; i<n; i++) {
    hionVector[i]->SetCutForSecondaryPhotons(cut);
    hionVector[i]->SetFluorescence(true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test17PhysicsList::SetCutForAugerElectrons(G4double cut)
{
  size_t n = hionVector.size();
  G4cout << " Cut for Auger electrons = " << cut/keV << " keV" << G4endl;

  for (size_t i=0; i<n; i++) {
    hionVector[i]->SetCutForAugerElectrons(cut);
    hionVector[i]->SetFluorescence(true);
    hionVector[i]->ActivateAugerElectronProduction(true);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....





