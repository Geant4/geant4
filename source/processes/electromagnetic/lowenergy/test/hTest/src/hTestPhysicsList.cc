// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//
// Always place it in front!
//
#include "G4Timer.hh"

#include "hTestPhysicsList.hh"
#include "hTestDetectorConstruction.hh"
#include "hTestPhysicsListMessenger.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4IonC12.hh"
#include "G4IonAr40.hh"
#include "G4IonFe56.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4EnergyLossTables.hh"
#include "G4ios.hh"
#include "g4std/iomanip"                

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPhysicsList::hTestPhysicsList(hTestDetectorConstruction* p)
:G4VUserPhysicsList()
{
  pDet = p;
 
  defaultCutValue = 0.1*cm;
  cutForGamma     = 1.0*cm;
  cutForElectron  = defaultCutValue;
  cutForProton    = defaultCutValue;
  
  //  MaxChargedStep = DBL_MAX; 
  MaxChargedStep = 0.1*cm; 
  
  SetVerboseLevel(2);
  physicsListMessenger = new hTestPhysicsListMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

hTestPhysicsList::~hTestPhysicsList()
{
  delete physicsListMessenger; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructBosons();
  ConstructLeptons();
  ConstructBarions();
  ConstructMesons();
  ConstructIons();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::ConstructBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  
  // gamma
  G4Gamma::GammaDefinition();
}
 //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::ConstructLeptons()
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

void hTestPhysicsList::ConstructMesons()
{
 //  mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::ConstructBarions()
{
//  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::ConstructIons()
{
//  Ions
  G4Deuteron::DeuteronDefinition();
  G4Alpha::AlphaDefinition();
  G4IonC12::IonC12Definition();
  G4IonAr40::IonAr40Definition();
  G4IonFe56::IonFe56Definition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructHad();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// gamma st
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

// gamma le
#include "G4LowEnergyRayleigh.hh"
#include "G4LowEnergyCompton.hh"   
#include "G4LowEnergyGammaConversion.hh"
#include "G4LowEnergyPhotoElectric.hh"


#include "G4MultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4LowEnergyIonisation.hh"
#include "G4LowEnergyBremsstrahlung.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCaptureAtRest.hh"

#include "G4hIonisation.hh"
#include "G4hLowEnergyIonisation.hh"

#include "hTestStepCut.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
      pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh);
      pmanager->AddDiscreteProcess(new G4LowEnergyPhotoElectric);
      pmanager->AddDiscreteProcess(new G4LowEnergyCompton);
      pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion);    

      //pmanager->AddDiscreteProcess(new G4PhotoElectricEffect);
      //pmanager->AddDiscreteProcess(new G4ComptonScattering);
      //pmanager->AddDiscreteProcess(new G4GammaConversion);    
      
    } else if (particleName == "e-") {
      //electron
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);

      //pmanager->AddProcess(new G4eIonisation,   -1, 2,2);
      //pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);       

       pmanager->AddProcess(new G4LowEnergyIonisation,  -1, 2,2);
       pmanager->AddProcess(new G4LowEnergyBremsstrahlung, -1,-1,3);   

      hTestStepCut* theeminusStepCut = new hTestStepCut();
      theeminusStepCut->SetMaxStep(MaxChargedStep);  
      pmanager->AddProcess(theeminusStepCut,         -1,-1,4);

    } else if (particleName == "e+") {
      //positron      
      pmanager->AddProcess(new G4MultipleScattering, -1, 1,1);
      //  pmanager->AddProcess(new G4LowEnergyIonisation, -1, 2,2);
      pmanager->AddProcess(new G4eIonisation,        -1, 2,2);
      pmanager->AddProcess(new G4eBremsstrahlung,    -1,-1,3);
      pmanager->AddProcess(new G4eplusAnnihilation,   0,-1,4);
                  
      hTestStepCut* theeplusStepCut = new hTestStepCut();
      theeplusStepCut->SetMaxStep(MaxChargedStep) ;          
      pmanager->AddProcess(theeplusStepCut,          -1,-1,5);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
     //muon  
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

      G4cout << "Hadronic processes for " << particleName << G4endl; 

      // Standard ionisation
      //            G4hIonisation* hIon = new G4hIonisation() ;

       G4hLowEnergyIonisation* hIon = new G4hLowEnergyIonisation() ;
      // hIon->SetNuclearStoppingOff() ;
      // hIon->SetNuclearStoppingOn() ;

      // hIon->SetStoppingPowerTableName("Ziegler1977He") ;
      // hIon->SetStoppingPowerTableName("Ziegler1977H") ;
      // hIon->SetStoppingPowerTableName("ICRU_R49p") ;
      // hIon->SetStoppingPowerTableName("ICRU_R49He") ;
      // hIon->SetStoppingPowerTableName("ICRU_R49PowersHe") ;

      // hIon->SetAntiProtonStoppingOn(); 

      pmanager->AddProcess(hIon,-1,2,2);
      
      hTestStepCut* thehadronStepCut = new hTestStepCut();
      thehadronStepCut->SetMaxStep(MaxChargedStep);          		       
      pmanager->AddProcess( thehadronStepCut,       -1,-1,3);
   
    } else if (   particleName == "alpha"  
               || particleName == "deuteron"  
               || particleName == "IonC12"  
               || particleName == "IonAr40"  
               || particleName == "IonFe56"  
              )
    {
      pmanager->AddProcess(new G4MultipleScattering,-1,1,1);

      G4cout << "Ionic processes for " << particleName << G4endl; 

      // Standard ionisation
      //   G4hIonisation* iIon = new G4hIonisation() ;

      // Standard ionisation with low energy extantion
      G4hLowEnergyIonisation* iIon = new G4hLowEnergyIonisation() ;

      // iIon->SetNuclearStoppingOff() ;
      // iIon->SetNuclearStoppingOn() ;
      // iIon->SetIonDefinition(particle) ;

      // iIon->SetStoppingPowerTableName("Ziegler1977He") ;
      // iIon->SetStoppingPowerTableName("Ziegler1977H") ;
      // iIon->SetStoppingPowerTableName("ICRU_R49p") ;
      // iIon->SetStoppingPowerTableName("ICRU_R49He") ;
      // iIon->SetStoppingPowerTableName("ICRU_R49PowersHe") ;

      pmanager->AddProcess(iIon,-1,2,2);
      
      hTestStepCut* theIonStepCut = new hTestStepCut();
      theIonStepCut->SetMaxStep(MaxChargedStep);          		       
      pmanager->AddProcess( theIonStepCut,       -1,-1,3);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// Hadron Processes

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4LambdaInelasticProcess.hh"
#include "G4AntiLambdaInelasticProcess.hh"
#include "G4SigmaPlusInelasticProcess.hh"
#include "G4SigmaMinusInelasticProcess.hh"
#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4XiZeroInelasticProcess.hh"
#include "G4XiMinusInelasticProcess.hh"
#include "G4AntiXiZeroInelasticProcess.hh"
#include "G4AntiXiMinusInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4OmegaMinusInelasticProcess.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"

// Low-energy Models

#include "G4LElastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"

#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::ConstructHad()
//
// Makes discrete physics processes for the hadrons, at present limited
// to those particles with GHEISHA interactions (INTRC > 0).
// The processes are: Elastic scattering, Inelastic scattering,
// Fission (for neutron only), and Capture (neutron).
//
// F.W.Jones  06-JUL-1998
//
{
   G4HadronElasticProcess* theElasticProcess = 
                                    new G4HadronElasticProcess;
   G4LElastic* theElasticModel = new G4LElastic;
   theElasticProcess->RegisterMe(theElasticModel);
   G4HadronElasticProcess* theElasticProcess1 = 
                                    new G4HadronElasticProcess;
   theParticleIterator->reset();
   while ((*theParticleIterator)()) {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "pi+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4PionPlusInelasticProcess* theInelasticProcess = 
                                new G4PionPlusInelasticProcess("inelastic");
         G4LEPionPlusInelastic* theInelasticModel = 
                                new G4LEPionPlusInelastic;
         theInelasticProcess->RegisterMe(theInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "pi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4PionMinusInelasticProcess* theInelasticProcess = 
                                new G4PionMinusInelasticProcess("inelastic");
         G4LEPionMinusInelastic* theInelasticModel = 
                                new G4LEPionMinusInelastic;
         theInelasticProcess->RegisterMe(theInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonPlusInelasticProcess* theInelasticProcess = 
                                  new G4KaonPlusInelasticProcess("inelastic");
         G4LEKaonPlusInelastic* theInelasticModel = new G4LEKaonPlusInelastic;
         theInelasticProcess->RegisterMe(theInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon0S") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonZeroSInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroSInelasticProcess("inelastic");
         G4LEKaonZeroSInelastic* theInelasticModel = 
                             new G4LEKaonZeroSInelastic;
         theInelasticProcess->RegisterMe(theInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon0L") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonZeroLInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroLInelasticProcess("inelastic");
         G4LEKaonZeroLInelastic* theInelasticModel = 
                             new G4LEKaonZeroLInelastic;
         theInelasticProcess->RegisterMe(theInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonMinusInelasticProcess* theInelasticProcess = 
                                 new G4KaonMinusInelasticProcess("inelastic");
         G4LEKaonMinusInelastic* theInelasticModel = 
                                 new G4LEKaonMinusInelastic;
         theInelasticProcess->RegisterMe(theInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4ProtonInelasticProcess* theInelasticProcess = 
                                    new G4ProtonInelasticProcess("inelastic");
         G4LEProtonInelastic* theInelasticModel = new G4LEProtonInelastic;
         theInelasticProcess->RegisterMe(theInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiProtonInelasticProcess* theInelasticProcess = 
                                new G4AntiProtonInelasticProcess("inelastic");
         G4LEAntiProtonInelastic* theInelasticModel = 
                                new G4LEAntiProtonInelastic;
         theInelasticProcess->RegisterMe(theInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "neutron") {
         
          // elastic scattering
         pmanager->AddDiscreteProcess(theElasticProcess);
         
          // inelastic scattering
         G4NeutronInelasticProcess* theInelasticProcess = 
                                    new G4NeutronInelasticProcess("inelastic");
         G4LENeutronInelastic* theInelasticModel = new G4LENeutronInelastic;
         theInelasticProcess->RegisterMe(theInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
         
	 /*
          // fission
         G4HadronFissionProcess* theFissionProcess =
                                    new G4HadronFissionProcess;
         G4LFission* theFissionModel = new G4LFission;
         theFissionProcess->RegisterMe(theFissionModel);
         pmanager->AddDiscreteProcess(theFissionProcess);
	 */         

         // capture
         G4HadronCaptureProcess* theCaptureProcess =
                                    new G4HadronCaptureProcess;
	 G4LCapture* theCaptureModel = new G4LCapture;
         theCaptureProcess->RegisterMe(theCaptureModel);
         pmanager->AddDiscreteProcess(theCaptureProcess);
      }
   }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4Decay.hh"

void hTestPhysicsList::ConstructGeneral()
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

void hTestPhysicsList::SetCuts()
{
  // G4Timer theTimer ;
  // theTimer.Start() ;

  //special for low energy physics
  //
  G4double lowlimit=1*eV;  
  G4Gamma   ::SetEnergyRange(lowlimit,100*GeV);
  G4Electron::SetEnergyRange(lowlimit,100*GeV);
  G4Positron::SetEnergyRange(lowlimit,100*GeV);
   
  if (verboseLevel >0){
    G4cout << "hTestPhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(MaxChargedStep,"Length") << G4endl;
  }  
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 

      G4cout << "Set cuts for all particles! " << G4endl; 

   SetCutValue(cutForGamma,"gamma");

   SetCutValue(cutForElectron,"e-");
   SetCutValue(cutForElectron,"e+");

   SetCutValue(cutForProton,"mu-");
   SetCutValue(cutForProton,"mu+");

  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton 

  SetCutValue(cutForProton, "proton");

  SetCutValue(cutForProton, "anti_proton");

  SetCutValueForOthers(cutForProton);

  //  G4cout << "Dump the table!" << G4endl;
              
  if (verboseLevel>0) DumpCutValuesTable();

  //  theTimer.Stop();
  //  G4cout.precision(6);
  //  G4cout << G4endl ;
  //  G4cout << "total time(SetCuts)=" << theTimer.GetUserElapsed() << " s " <<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::SetGammaCut(G4double val)
{
  ResetCuts();
  cutForGamma = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::SetElectronCut(G4double val)
{
  ResetCuts();
  cutForElectron = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::SetProtonCut(G4double val)
{
  ResetCuts();
  cutForProton = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::SetCutsByEnergy(G4double val)
{
  G4ParticleTable* thehTestParticleTable =  G4ParticleTable::GetParticleTable();
  G4Material* currMat = pDet->GetAbsorberMaterial();

  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma
  G4ParticleDefinition* part;
  G4double cut;

  part = thehTestParticleTable->FindParticle("e-");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "e-");

  part = thehTestParticleTable->FindParticle("e+");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "e+");

  // set cut values for proton and anti_proton before all other hadrons
  // because some processes for hadrons need cut values for proton/anti_proton
  part = thehTestParticleTable->FindParticle("proton");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "proton");

  part = thehTestParticleTable->FindParticle("anti_proton");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  SetCutValue(cut, "anti_proton");

  SetCutValueForOthers(cut);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::GetRange(G4double val)
{
  G4ParticleTable* thehTestParticleTable =  G4ParticleTable::GetParticleTable();
  G4Material* currMat = pDet->GetAbsorberMaterial();

  G4ParticleDefinition* part;
  G4double cut;
  part = thehTestParticleTable->FindParticle("e-");
  cut = G4EnergyLossTables::GetRange(part,val,currMat);
  G4cout << "material : " << currMat->GetName() << G4endl;
  G4cout << "particle : " << part->GetParticleName() << G4endl;
  G4cout << "energy   : " << val / keV << " (keV)" << G4endl;
  G4cout << "range    : " << cut / mm << " (mm)" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void hTestPhysicsList::SetMaxStep(G4double step)
{
  MaxChargedStep = step ;
  G4cout << " MaxChargedStep=" << MaxChargedStep << G4endl;
  G4cout << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....







