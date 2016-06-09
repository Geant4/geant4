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
#include "ExGflashPhysicsList.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ios.hh"
#include <iomanip>   

#include "G4FastSimulationManagerProcess.hh"

using namespace std;

ExGflashPhysicsList::ExGflashPhysicsList():  G4VUserPhysicsList()
{
  SetVerboseLevel(0);
}

ExGflashPhysicsList::~ExGflashPhysicsList()
{
}

void ExGflashPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

	std::cout<<"start construct particle"<<std::endl;
  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
  ConstructIons();
	std::cout<<"end construct particle"<<std::endl;
}

void ExGflashPhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

#include "G4LeptonConstructor.hh"
void ExGflashPhysicsList::ConstructLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

#include "G4MesonConstructor.hh"
void ExGflashPhysicsList::ConstructMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

#include "G4BaryonConstructor.hh"
void ExGflashPhysicsList::ConstructBaryons()
{
  //  Construct all barions
  G4BaryonConstructor  pConstructor;
  pConstructor.ConstructParticle();  
}

#include "G4IonConstructor.hh"
void ExGflashPhysicsList::ConstructIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void ExGflashPhysicsList::ConstructProcess()
{
  //	std::cout<<"1111"<<std::endl;
  AddTransportation();
  // 	std::cout<<"2222"<<std::endl;
  AddParameterisation();
  std::cout<<"AddParameterisation"<<std::endl;
 
  ConstructEM();
  std::cout<<"ConstructEM"<<std::endl;
  ConstructGeneral();
  //	std::cout<<"5555"<<std::endl;
}


void ExGflashPhysicsList::AddTransportation()
{
  G4VUserPhysicsList::AddTransportation();
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4UserLimits.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"
void ExGflashPhysicsList::ConstructEM()
{

   G4cout<<"Physics List constructor"<<G4endl;
   SetCuts();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
      // Construct processes for gamma
      G4VProcess* theGammaConversion = new G4GammaConversion();
      G4VProcess* theComptonScattering = new G4ComptonScattering();
      G4VProcess* thePhotoElectricEffect = new G4PhotoElectricEffect();     
      // G4VProcess* thegammacut = new G4UserLimits();
      //      thegammacut->SetUserMinEkine(1.0*MeV);

      pmanager->AddDiscreteProcess(theGammaConversion);
      pmanager->AddDiscreteProcess(theComptonScattering);      
      pmanager->AddDiscreteProcess(thePhotoElectricEffect);
      //  pmanager->AddProcess(thegammacut);
      //   G4cout <<"theGammaConversion" << theGammaConversion <<endl;
      //G4cout <<"theComptonScattering" << theComptonScattering <<endl;
      //G4cout <<"thePhotoElectricEffect" << thePhotoElectricEffect <<endl;

    } else if (particleName == "e-") {
    //electron
      // Construct processes for electron
      G4VProcess* theeminusMultipleScattering = new G4eMultipleScattering();
      G4VProcess* theeminusIonisation = new G4eIonisation();
      G4VProcess* theeminusBremsstrahlung = new G4eBremsstrahlung();
      //      G4VProcess* theeminuscut = new G4UserLimits();
      //  theeminuscut->SetUserMinEkine(1.0*MeV);
      // add processes
      pmanager->AddProcess(theeminusMultipleScattering);
      pmanager->AddProcess(theeminusIonisation);
      pmanager->AddProcess(theeminusBremsstrahlung);      
     
      //  pmanager->AddProcess( theeminuscut);


  
 // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxAlongStep,  1);
      pmanager->SetProcessOrdering(theeminusIonisation, idxAlongStep,  2);
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(theeminusIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung, idxPostStep, 3);

      //G4cout <<"theeminusMultipleScattering" << theeminusMultipleScattering <<endl;
      //G4cout <<"theeminusIonisation" << theeminusIonisation <<endl;
      //G4cout <<"theeminusBremsstrahlung" << theeminusBremsstrahlung <<endl;

    } else if (particleName == "e+") {
    //positron
      // Construct processes for positron
      G4VProcess* theeplusMultipleScattering = new G4eMultipleScattering();
      G4VProcess* theeplusIonisation = new G4eIonisation();
      G4VProcess* theeplusBremsstrahlung = new G4eBremsstrahlung();
      G4VProcess* theeplusAnnihilation = new G4eplusAnnihilation();
      // add processes
      pmanager->AddProcess(theeplusMultipleScattering);
      pmanager->AddProcess(theeplusIonisation);
      pmanager->AddProcess(theeplusBremsstrahlung);
      pmanager->AddProcess(theeplusAnnihilation);
      // set ordering for AtRestDoIt
      pmanager->SetProcessOrderingToFirst(theeplusAnnihilation, idxAtRest);
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxAlongStep,  1);
      pmanager->SetProcessOrdering(theeplusIonisation, idxAlongStep,  2);
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeplusMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(theeplusIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(theeplusBremsstrahlung, idxPostStep, 3);
      pmanager->SetProcessOrdering(theeplusAnnihilation, idxPostStep, 4);
  
      //G4cout <<"theeplusMultipleScattering" << theeplusMultipleScattering <<endl;
      //G4cout <<"theeplusIonisation" << theeplusIonisation <<endl;
      //G4cout <<"theeplusBremsstrahlung" << theeplusBremsstrahlung <<endl;

    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
    //muon  
     // Construct processes for muon+
     G4VProcess* aMultipleScattering = new G4MuMultipleScattering();
     G4VProcess* aBremsstrahlung = new G4MuBremsstrahlung();
     G4VProcess* aPairProduction = new G4MuPairProduction();
     G4VProcess* anIonisation = new G4MuIonisation();
      // add processes
     pmanager->AddProcess(anIonisation);
     pmanager->AddProcess(aMultipleScattering);
     pmanager->AddProcess(aBremsstrahlung);
     pmanager->AddProcess(aPairProduction);
     // set ordering for AlongStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,  1);
     pmanager->SetProcessOrdering(anIonisation, idxAlongStep,  2);
     // set ordering for PostStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);
     pmanager->SetProcessOrdering(anIonisation, idxPostStep, 2);
     pmanager->SetProcessOrdering(aBremsstrahlung, idxPostStep, 3);
     pmanager->SetProcessOrdering(aPairProduction, idxPostStep, 4);
     
   } else if ((!particle->IsShortLived()) &&
	      (particle->GetPDGCharge() != 0.0) && 
	      (particle->GetParticleName() != "chargedgeantino")) {
     // all others charged particles except geantino
     G4VProcess* aMultipleScattering = new G4hMultipleScattering();
     G4VProcess* anIonisation = new G4hIonisation();
     // add processes
     pmanager->AddProcess(anIonisation);
     pmanager->AddProcess(aMultipleScattering);
     // set ordering for AlongStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,  1);
     pmanager->SetProcessOrdering(anIonisation, idxAlongStep,  2);
     // set ordering for PostStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);
     pmanager->SetProcessOrdering(anIonisation, idxPostStep, 2);
    }
  }
}


#include "G4Decay.hh"
void ExGflashPhysicsList::ConstructGeneral()
{
  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  //G4cout << "decay" <<theDecayProcess<<endl;
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

//  WARNING: This methode is mandatory if U want to use GFLASH
void ExGflashPhysicsList::AddParameterisation()
{
  G4FastSimulationManagerProcess* 
    theFastSimulationManagerProcess = 
      new G4FastSimulationManagerProcess();
  G4cout << "FastSimulationManagerProcess" <<G4endl;
  theParticleIterator->reset();
  //std::cout<<"---"<<std::endl;
  while( (*theParticleIterator)() ){
    //std::cout<<"+++"<<std::endl;
  
    G4ParticleDefinition* particle = theParticleIterator->value();
    // std::cout<<"--- particle "<<particle->GetParticleName()<<std::endl;
    G4ProcessManager* pmanager = particle->GetProcessManager();
    // The fast simulation process becomes a discrete process only since 9.0:
    pmanager->AddDiscreteProcess(theFastSimulationManagerProcess);
  }
}

void ExGflashPhysicsList::SetCuts()
{
  if (verboseLevel >1){
    G4cout << "ExGflashPhysicsList::SetCuts:";
  }  
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  DumpCutValuesTable();
   SetCutsWithDefault();
   //  SetCutValue(100*mm, "gamma");
//   SetCutValue(0*mm, "e-");
//   SetCutValue(0*mm, "e+");

//   SetCutValue(62*mm, "gamma");
//   SetCutValue(0.73*mm, "e-");
//   SetCutValue(0.78*mm, "e+");



 
 
  DumpCutValuesTable();
// set cuts for region crystals with default Cuts
  G4Region* region = G4RegionStore::GetInstance()->GetRegion("crystals");
  region->SetProductionCuts(
          G4ProductionCutsTable::GetProductionCutsTable()->GetDefaultProductionCuts());
}









