// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T08PhysicsList.cc,v 1.2 1999-04-17 07:24:03 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// ------------------------------------------------------------
//	GEANT 4 class header file 
//
//      This is a version for maximum particle set
// ------------------------------------------------------------

#include "globals.hh"
#include "T08PhysicsList.hh"
#include "T08PhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip.h>                


T08PhysicsList::T08PhysicsList():  G4VUserPhysicsList(),
 thePhotoElectricEffect(NULL),theComptonScattering(NULL),
 theGammaConversion(NULL),
 theeminusMultipleScattering(NULL),theeminusIonisation(NULL),
 theeminusBremsstrahlung(NULL),
 theeplusMultipleScattering(NULL),theeplusIonisation(NULL),
 theeplusBremsstrahlung(NULL),
 theeplusAnnihilation(NULL)
{
  defaultCutValue = 2.0*mm;
  SetVerboseLevel(1);
  physicsListMessenger = new T08PhysicsListMessenger(this);
}

T08PhysicsList::~T08PhysicsList()
{
}

void T08PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBarions();

}

void T08PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}

void T08PhysicsList::ConstructLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

void T08PhysicsList::ConstructMesons()
{
 //  mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4RhoZero::RhoZeroDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
}

void T08PhysicsList::ConstructBarions()
{
//  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}


void T08PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

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

#include "G4hIonisation.hh"
///#include "G4UserSpecialCuts.hh"

void T08PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
      // Construct processes for gamma
      thePhotoElectricEffect = new G4PhotoElectricEffect();      
      theComptonScattering = new G4ComptonScattering();
      theGammaConversion = new G4GammaConversion();
      
      pmanager->AddDiscreteProcess(thePhotoElectricEffect);
      pmanager->AddDiscreteProcess(theComptonScattering);
      pmanager->AddDiscreteProcess(theGammaConversion);
      
    } else if (particleName == "e-") {
    //electron
      // Construct processes for electron
      theeminusMultipleScattering = new G4MultipleScattering();
      theeminusIonisation = new G4eIonisation();
      theeminusBremsstrahlung = new G4eBremsstrahlung();
      // add processes
      pmanager->AddProcess(theeminusMultipleScattering);
      pmanager->AddProcess(theeminusIonisation);
      pmanager->AddProcess(theeminusBremsstrahlung);      
      // set ordering for AlongStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxAlongStep,  1);
      pmanager->SetProcessOrdering(theeminusIonisation, idxAlongStep,  2);
      // set ordering for PostStepDoIt
      pmanager->SetProcessOrdering(theeminusMultipleScattering, idxPostStep, 1);
      pmanager->SetProcessOrdering(theeminusIonisation, idxPostStep, 2);
      pmanager->SetProcessOrdering(theeminusBremsstrahlung, idxPostStep, 3);

    } else if (particleName == "e+") {
    //positron
      // Construct processes for positron
      theeplusMultipleScattering = new G4MultipleScattering();
      theeplusIonisation = new G4eIonisation();
      theeplusBremsstrahlung = new G4eBremsstrahlung();
      theeplusAnnihilation = new G4eplusAnnihilation();
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
  
    } 
#ifdef ENABLE_MUON_PROCESSES
    else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
    //muon  
     // Construct processes for muon+
     G4VProcess* aMultipleScattering = new G4MultipleScattering();
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

    }
#endif
      else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
     // all others charged particles except geantino
     G4VProcess* aMultipleScattering = new G4MultipleScattering();
     G4VProcess* anIonisation = new G4hIonisation();
     //>>G4VProcess*  theUserSpecialCuts = new G4UserSpecialCuts();
     // add processes
     pmanager->AddProcess(anIonisation);
     pmanager->AddProcess(aMultipleScattering);
     /// pmanager->AddProcess(theUserSpecialCuts);
     // set ordering for AtRestDoIt
     //>>pmanager->SetProcessOrderingToFirst(theUserSpecialCuts, idxAtRest);
     // set ordering for AlongStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,  1);
     pmanager->SetProcessOrdering(anIonisation, idxAlongStep,  2);
     // set ordering for PostStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);
     pmanager->SetProcessOrdering(anIonisation, idxPostStep, 2);
     //>>pmanager->SetProcessOrdering(theUserSpecialCuts, idxPostStep, 3);
    }
  }
}

#include "G4Decay.hh"
void T08PhysicsList::ConstructGeneral()
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

void T08PhysicsList::SetCuts()
{
  G4double cut = defaultCutValue;
  if (verboseLevel >0){
    G4cout << "T08PhysicsList::SetCuts:";
    G4cout << "CutLength : " << cut/mm << " (mm)" << endl;
  }  

 //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}
 


void T08PhysicsList::SetStatusEmProcess()
{
// Process definitions

   const
   G4String process1("phot"),process2("Comp"),process3("conv"),
            process4("muls"),process5("ioni"),process6("brem"),
            process7("anni");
   static        
   G4String lstate1("on"),lstate2("on"),lstate3("on"),
            lstate4("on"),lstate5("on"),lstate6("on"),lstate7("on");
             
   G4String process,newstate;
   G4bool active(true);

   G4ProcessManager* pGammaManager = G4Gamma::GammaDefinition()       -> GetProcessManager();
   G4ProcessManager* pElectManager = G4Electron::ElectronDefinition() -> GetProcessManager();
   G4ProcessManager* pPositManager = G4Positron::PositronDefinition() -> GetProcessManager();
//
//
// Loop over processes
//
 do {
     // diplay the status of the processes
        G4cout << " ****** Processes status ****** " << endl;
        G4cout << "   " << process1 << " " << process2 << " " << process3
             <<  " "  << process4 << " " << process5 << " " << process6 
             <<  " "  << process7 << endl;
        G4cout << "   "  << lstate1 << "   " << lstate2 << "   " << lstate3
             << "   " << lstate4 << "   " << lstate5 << "   " << lstate6 
             << "   " << lstate7 << endl;

     // update the status of the processes
        G4cout << "  enter -> a process name, on/off or: ok on : " << flush;
        cin >> process >> newstate ;
        if (newstate=="on") active = true;
         else active = false;

        if (process == "phot")
           { pGammaManager -> SetProcessActivation (thePhotoElectricEffect, active);
             lstate1 = newstate;
           }
        else if (process == "Comp")
           { pGammaManager -> SetProcessActivation (theComptonScattering, active);
             lstate2 = newstate;
           }
        else if (process == "conv")
           { pGammaManager -> SetProcessActivation (theGammaConversion, active);
             lstate3 = newstate;
           }
        else if (process == "muls")
           { pElectManager -> SetProcessActivation (theeminusMultipleScattering, active);
             pPositManager -> SetProcessActivation (theeplusMultipleScattering, active);
             lstate4 = newstate;
           }               
        else if (process == "ioni")
           { pElectManager -> SetProcessActivation (theeminusIonisation, active);
             pPositManager -> SetProcessActivation (theeplusIonisation, active);
             lstate5 = newstate;
           }
        else if (process == "brem")
           { pElectManager -> SetProcessActivation (theeminusBremsstrahlung, active);
             pPositManager -> SetProcessActivation (theeplusBremsstrahlung, active);
             lstate6 = newstate;
           }
        else if (process == "anni")
           { pPositManager -> SetProcessActivation (theeplusAnnihilation, active);
             lstate7 = newstate;
           }
  } while (process!="ok") ;
  
 G4cout << "  ---> Done" << endl;  
}




