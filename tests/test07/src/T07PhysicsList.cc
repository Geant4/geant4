// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: T07PhysicsList.cc,v 1.5 2000-02-25 13:32:17 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "T07PhysicsList.hh"
#include "T07PhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4ios.hh"
#include "g4std/iomanip"                

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

T07PhysicsList::T07PhysicsList():  G4VUserPhysicsList(),
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
  physicsListMessenger = new T07PhysicsListMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

T07PhysicsList::~T07PhysicsList()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructMesons();
  ConstructBaryons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07PhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();

  // gamma
  G4Gamma::GammaDefinition();

  // optical photon
  G4OpticalPhoton::OpticalPhotonDefinition();
}
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07PhysicsList::ConstructLeptons()
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07PhysicsList::ConstructMesons()
{
 //  mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07PhysicsList::ConstructBaryons()
{
//  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07PhysicsList::ConstructEM()
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

      pmanager->AddProcess(theeminusMultipleScattering,-1,1,1);
      pmanager->AddProcess(theeminusIonisation,-1,2,2);
      pmanager->AddProcess(theeminusBremsstrahlung,-1,-1,3);      

    } else if (particleName == "e+") {
    //positron
      // Construct processes for positron
      theeplusMultipleScattering = new G4MultipleScattering();
      theeplusIonisation = new G4eIonisation();
      theeplusBremsstrahlung = new G4eBremsstrahlung();
      theeplusAnnihilation = new G4eplusAnnihilation();
      
      pmanager->AddProcess(theeplusMultipleScattering,-1,1,1);
      pmanager->AddProcess(theeplusIonisation,-1,2,2);
      pmanager->AddProcess(theeplusBremsstrahlung,-1,-1,3);
      pmanager->AddProcess(theeplusAnnihilation,0,-1,4);      
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
    //muon  
     // Construct processes for muon+
     pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
     pmanager->AddProcess(new G4MuIonisation(),-1,2,2);
     pmanager->AddProcess(new G4MuBremsstrahlung(),-1,-1,3);
     pmanager->AddProcess(new G4MuPairProduction(),-1,-1,4);       
     
    } else if ((!particle->IsShortLived()) &&
	       (particle->GetPDGCharge() != 0.0) && 
	       (particle->GetParticleName() != "chargedgeantino")) {
     // all others charged particles except geantino
       pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
       pmanager->AddProcess(new G4hIonisation(),-1,2,2);      
    }
  }
}


#include "G4Decay.hh"

void T07PhysicsList::ConstructGeneral()
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

void T07PhysicsList::SetCuts()
{
  G4double cut = defaultCutValue;
  if (verboseLevel >0){
    G4cout << "T07PhysicsList::SetCuts:";
    G4cout << "CutLength : " << cut/mm << " (mm)" << G4endl;
  }  

 //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}
 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void T07PhysicsList::SetStatusEmProcess()
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
        G4cout << " ****** Processes status ****** " << G4endl;
        G4cout << "   " << process1 << " " << process2 << " " << process3
             <<  " "  << process4 << " " << process5 << " " << process6 
             <<  " "  << process7 << G4endl;
        G4cout << "   "  << lstate1 << "   " << lstate2 << "   " << lstate3
             << "   " << lstate4 << "   " << lstate5 << "   " << lstate6 
             << "   " << lstate7 << G4endl;

     // update the status of the processes
        G4cout << "  enter -> a process name, on/off or: ok on : " << G4std::flush;
        G4cin >> process >> newstate ;
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
  
 G4cout << "  ---> Done" << G4endl;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....




