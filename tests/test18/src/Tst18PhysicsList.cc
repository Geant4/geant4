// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Tst18PhysicsList.cc,v 1.3 2000-06-06 12:15:08 flei Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "Tst18PhysicsList.hh"

#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4ios.hh"
#include <iomanip.h>   

#include "G4IonTable.hh"
#include "G4Ions.hh"

#include "G4RadioactiveDecay.hh"


//#include "G4FastSimulationManagerProcess.hh"


Tst18PhysicsList::Tst18PhysicsList():  G4VUserPhysicsList()
{
  // default cut value  (1.0mm) 
  defaultCutValue = 1.*m;

  SetVerboseLevel(1);
}

Tst18PhysicsList::~Tst18PhysicsList()
{
}

void Tst18PhysicsList::ConstructParticle()
{

  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  // create all particles
  ConstructAllBosons();
  ConstructAllLeptons();
  ConstructAllMesons();
  ConstructAllBaryons();
  ConstructAllIons();
  ConstructAllShortLiveds();

  //  G4GenericIon::GenericIonDefinition();
}

void Tst18PhysicsList::ConstructProcess()
{
  AddTransportation();

  ConstructEM();
  ConstructHad();
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
#include "G4ionIonisation.hh"

void Tst18PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
      // Construct processes for gamma
      //      pmanager->AddDiscreteProcess(new G4GammaConversion());
      //pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      //pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
    //electron
      // Construct processes for electron
      G4VProcess* theeminusMultipleScattering = new G4MultipleScattering();
      G4VProcess* theeminusIonisation = new G4eIonisation();
      G4VProcess* theeminusBremsstrahlung = new G4eBremsstrahlung();
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
      G4VProcess* theeplusMultipleScattering = new G4MultipleScattering();
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
  
    } else if( particleName == "mu+" || 
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
     
     /*    } else if( particleName == "GenericIon" ) {
     G4VProcess* aionIonization = new G4ionIonisation;
     G4VProcess* aMultipleScattering = new G4MultipleScattering();
     pmanager->AddProcess(aionIonization);
     pmanager->AddProcess(aMultipleScattering);
     // set ordering for AlongStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxAlongStep,  1);
     pmanager->SetProcessOrdering(aionIonization, idxAlongStep,  2);
     // set ordering for PostStepDoIt
     pmanager->SetProcessOrdering(aMultipleScattering, idxPostStep, 1);
     pmanager->SetProcessOrdering(aionIonization, idxPostStep, 2); 
  */
   } else if ((!particle->IsShortLived()) &&
	      (particle->GetPDGCharge() != 0.0) && 
	      (particle->GetParticleName() != "chargedgeantino")) {
     // all others charged particles except geantino

      pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
      G4hIonisation*    hIon;
#if 0 
      if ( 0 ) { 
	G4hEnergyLoss* hIonLowE;

	// hIon=  new G4hLowEnergyIonisation() ;
	// Definition of nuclear stopping
        if( 1 ) { 
	  hIonLowE->SetNuclearStoppingOff() ;
	  // Definition of the proton's table  
	  hIonLowE->SetStoppingPowerTableName("ICRU_R49p") ;

        } else { // to use Ziegler: 
	  hIonLowE->SetNuclearStoppingOn() ;
	  // Definition of the proton's table  
	  hIonLowE->SetStoppingPowerTableName("Ziegler1977H") ;
	}
        hIon= hIonLowE;
      }else{
	hIon= new G4hIonisation();
      }
#endif
      hIon= new G4hIonisation();
      pmanager->AddProcess(hIon,  -1,2,2);      

      G4double demax = 0.05;  // try to lose at most 5% of the energy in 
                              //    a single step (in limit of large energies)
      G4double stmin = 0.01 * m;  // length of the final step 
      hIon->SetStepFunction( demax, stmin );
    }
  }
}


// Hadron Processes

#include "G4HadronElasticProcess.hh"

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

#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LEAntiNeutronInelastic.hh"
#include "G4LELambdaInelastic.hh"
#include "G4LEAntiLambdaInelastic.hh"
#include "G4LESigmaPlusInelastic.hh"
#include "G4LESigmaMinusInelastic.hh"
#include "G4LEAntiSigmaPlusInelastic.hh"
#include "G4LEAntiSigmaMinusInelastic.hh"
#include "G4LEXiZeroInelastic.hh"
#include "G4LEXiMinusInelastic.hh"
#include "G4LEAntiXiZeroInelastic.hh"
#include "G4LEAntiXiMinusInelastic.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4LEOmegaMinusInelastic.hh"
#include "G4LEAntiOmegaMinusInelastic.hh"

// High-energy Models

#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"
#include "G4HEProtonInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4HENeutronInelastic.hh"
#include "G4HEAntiNeutronInelastic.hh"
#include "G4HELambdaInelastic.hh"
#include "G4HEAntiLambdaInelastic.hh"
#include "G4HESigmaPlusInelastic.hh"
#include "G4HESigmaMinusInelastic.hh"
#include "G4HEAntiSigmaPlusInelastic.hh"
#include "G4HEAntiSigmaMinusInelastic.hh"
#include "G4HEXiZeroInelastic.hh"
#include "G4HEXiMinusInelastic.hh"
#include "G4HEAntiXiZeroInelastic.hh"
#include "G4HEAntiXiMinusInelastic.hh"
#include "G4HEOmegaMinusInelastic.hh"
#include "G4HEAntiOmegaMinusInelastic.hh"

// Stopping processes

#ifdef TRIUMF_STOP_PIMINUS
#include "G4PionMinusAbsorptionAtRest.hh"
#else
#include "G4PiMinusAbsorptionAtRest.hh"
#endif
#ifdef TRIUMF_STOP_KMINUS
#include "G4KaonMinusAbsorption.hh"
#else
#include "G4KaonMinusAbsorptionAtRest.hh"
#endif

//
// ConstructHad()
//
// Makes discrete physics processes for the hadrons, at present limited
// to those particles with GHEISHA interactions (INTRC > 0).
// The processes are: Elastic scattering and Inelastic scattering.
//
// F.W.Jones  09-JUL-1998
//

void Tst18PhysicsList::ConstructHad()
{
   G4HadronElasticProcess* theElasticProcess = 
                                    new G4HadronElasticProcess;
   G4LElastic* theElasticModel = new G4LElastic;
   theElasticProcess->RegisterMe(theElasticModel);

   theParticleIterator->reset();
   while ((*theParticleIterator)()) {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "pi+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4PionPlusInelasticProcess* theInelasticProcess = 
                                new G4PionPlusInelasticProcess("inelastic");
         G4LEPionPlusInelastic* theLEInelasticModel = 
                                new G4LEPionPlusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEPionPlusInelastic* theHEInelasticModel = 
                                new G4HEPionPlusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "pi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4PionMinusInelasticProcess* theInelasticProcess = 
                                new G4PionMinusInelasticProcess("inelastic");
         G4LEPionMinusInelastic* theLEInelasticModel = 
                                new G4LEPionMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEPionMinusInelastic* theHEInelasticModel = 
                                new G4HEPionMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
#ifdef TRIUMF_STOP_PIMINUS
         pmanager->AddRestProcess(new G4PionMinusAbsorptionAtRest, ordDefault);
#else
         G4String prcNam;
         pmanager->AddRestProcess(
           new G4PiMinusAbsorptionAtRest(
                prcNam="PiMinusAbsorptionAtRest"), ordDefault);
#endif
      }
      else if (particleName == "kaon+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonPlusInelasticProcess* theInelasticProcess = 
                                  new G4KaonPlusInelasticProcess("inelastic");
         G4LEKaonPlusInelastic* theLEInelasticModel = 
                                  new G4LEKaonPlusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEKaonPlusInelastic* theHEInelasticModel = 
                                  new G4HEKaonPlusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon0S") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonZeroSInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroSInelasticProcess("inelastic");
         G4LEKaonZeroSInelastic* theLEInelasticModel = 
                             new G4LEKaonZeroSInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEKaonZeroInelastic* theHEInelasticModel = 
                             new G4HEKaonZeroInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon0L") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonZeroLInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroLInelasticProcess("inelastic");
         G4LEKaonZeroLInelastic* theLEInelasticModel = 
                             new G4LEKaonZeroLInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEKaonZeroInelastic* theHEInelasticModel = 
                             new G4HEKaonZeroInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonMinusInelasticProcess* theInelasticProcess = 
                                 new G4KaonMinusInelasticProcess("inelastic");
         G4LEKaonMinusInelastic* theLEInelasticModel = 
                                 new G4LEKaonMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEKaonMinusInelastic* theHEInelasticModel = 
                                 new G4HEKaonMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
#ifdef TRIUMF_STOP_KMINUS
         pmanager->AddRestProcess(new G4KaonMinusAbsorption, ordDefault);
#else
         pmanager->AddRestProcess(new G4KaonMinusAbsorptionAtRest, ordDefault);
#endif
      }
      else if (particleName == "proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4ProtonInelasticProcess* theInelasticProcess = 
                                    new G4ProtonInelasticProcess("inelastic");
         G4LEProtonInelastic* theLEInelasticModel = new G4LEProtonInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEProtonInelastic* theHEInelasticModel = new G4HEProtonInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiProtonInelasticProcess* theInelasticProcess = 
                                new G4AntiProtonInelasticProcess("inelastic");
         G4LEAntiProtonInelastic* theLEInelasticModel = 
                                new G4LEAntiProtonInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiProtonInelastic* theHEInelasticModel = 
                                new G4HEAntiProtonInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "neutron") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4NeutronInelasticProcess* theInelasticProcess = 
                                    new G4NeutronInelasticProcess("inelastic");
         G4LENeutronInelastic* theLEInelasticModel = 
                                    new G4LENeutronInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HENeutronInelastic* theHEInelasticModel = 
                                    new G4HENeutronInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }  
      else if (particleName == "anti_neutron") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiNeutronInelasticProcess* theInelasticProcess = 
                               new G4AntiNeutronInelasticProcess("inelastic");
         G4LEAntiNeutronInelastic* theLEInelasticModel = 
                               new G4LEAntiNeutronInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiNeutronInelastic* theHEInelasticModel = 
                               new G4HEAntiNeutronInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "lambda") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4LambdaInelasticProcess* theInelasticProcess = 
                                    new G4LambdaInelasticProcess("inelastic");
         G4LELambdaInelastic* theLEInelasticModel = new G4LELambdaInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HELambdaInelastic* theHEInelasticModel = new G4HELambdaInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_lambda") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiLambdaInelasticProcess* theInelasticProcess = 
                                new G4AntiLambdaInelasticProcess("inelastic");
         G4LEAntiLambdaInelastic* theLEInelasticModel = 
                                new G4LEAntiLambdaInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiLambdaInelastic* theHEInelasticModel = 
                                new G4HEAntiLambdaInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "sigma+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4SigmaPlusInelasticProcess* theInelasticProcess = 
                                 new G4SigmaPlusInelasticProcess("inelastic");
         G4LESigmaPlusInelastic* theLEInelasticModel = 
                                 new G4LESigmaPlusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HESigmaPlusInelastic* theHEInelasticModel = 
                                 new G4HESigmaPlusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "sigma-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4SigmaMinusInelasticProcess* theInelasticProcess = 
                                 new G4SigmaMinusInelasticProcess("inelastic");
         G4LESigmaMinusInelastic* theLEInelasticModel = 
                                 new G4LESigmaMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HESigmaMinusInelastic* theHEInelasticModel = 
                                 new G4HESigmaMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_sigma+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiSigmaPlusInelasticProcess* theInelasticProcess = 
                             new G4AntiSigmaPlusInelasticProcess("inelastic");
         G4LEAntiSigmaPlusInelastic* theLEInelasticModel = 
                                 new G4LEAntiSigmaPlusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiSigmaPlusInelastic* theHEInelasticModel = 
                                 new G4HEAntiSigmaPlusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_sigma-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiSigmaMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiSigmaMinusInelasticProcess("inelastic");
         G4LEAntiSigmaMinusInelastic* theLEInelasticModel = 
                                 new G4LEAntiSigmaMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiSigmaMinusInelastic* theHEInelasticModel = 
                                 new G4HEAntiSigmaMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "xi0") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4XiZeroInelasticProcess* theInelasticProcess = 
                            new G4XiZeroInelasticProcess("inelastic");
         G4LEXiZeroInelastic* theLEInelasticModel = 
                                 new G4LEXiZeroInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEXiZeroInelastic* theHEInelasticModel = 
                                 new G4HEXiZeroInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "xi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4XiMinusInelasticProcess* theInelasticProcess = 
                            new G4XiMinusInelasticProcess("inelastic");
         G4LEXiMinusInelastic* theLEInelasticModel = 
                                 new G4LEXiMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEXiMinusInelastic* theHEInelasticModel = 
                                 new G4HEXiMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_xi0") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiXiZeroInelasticProcess* theInelasticProcess = 
                            new G4AntiXiZeroInelasticProcess("inelastic");
         G4LEAntiXiZeroInelastic* theLEInelasticModel = 
                                 new G4LEAntiXiZeroInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiXiZeroInelastic* theHEInelasticModel = 
                                 new G4HEAntiXiZeroInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_xi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiXiMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiXiMinusInelasticProcess("inelastic");
         G4LEAntiXiMinusInelastic* theLEInelasticModel = 
                                 new G4LEAntiXiMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiXiMinusInelastic* theHEInelasticModel = 
                                 new G4HEAntiXiMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "deuteron") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4DeuteronInelasticProcess* theInelasticProcess = 
                            new G4DeuteronInelasticProcess("inelastic");
         G4LEDeuteronInelastic* theLEInelasticModel = 
                                 new G4LEDeuteronInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "triton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4TritonInelasticProcess* theInelasticProcess = 
                            new G4TritonInelasticProcess("inelastic");
         G4LETritonInelastic* theLEInelasticModel = 
                                 new G4LETritonInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "alpha") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AlphaInelasticProcess* theInelasticProcess = 
                            new G4AlphaInelasticProcess("inelastic");
         G4LEAlphaInelastic* theLEInelasticModel = 
                                 new G4LEAlphaInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "omega-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4OmegaMinusInelasticProcess* theInelasticProcess = 
                            new G4OmegaMinusInelasticProcess("inelastic");
         G4LEOmegaMinusInelastic* theLEInelasticModel = 
                                 new G4LEOmegaMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEOmegaMinusInelastic* theHEInelasticModel = 
                                 new G4HEOmegaMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_omega-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiOmegaMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiOmegaMinusInelasticProcess("inelastic");
         G4LEAntiOmegaMinusInelastic* theLEInelasticModel = 
                                 new G4LEAntiOmegaMinusInelastic;
         theInelasticProcess->RegisterMe(theLEInelasticModel);
         G4HEAntiOmegaMinusInelastic* theHEInelasticModel = 
                                 new G4HEAntiOmegaMinusInelastic;
         theInelasticProcess->RegisterMe(theHEInelasticModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
   }
}


#include "G4Decay.hh"
void Tst18PhysicsList::ConstructGeneral()
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
  // Declare radioactive decay to the GenericIon in the IonTable.
  //
  const G4IonTable *theIonTable = G4ParticleTable::GetParticleTable()->GetIonTable();
  G4RadioactiveDecay *theRadioactiveDecay = new G4RadioactiveDecay();
  for (G4int i=0; i<theIonTable->Entries(); i++)
    {
      G4String particleName = theIonTable->GetParticle(i)->GetParticleName();
      if (particleName == "GenericIon")
	{
	  G4ProcessManager* pmanager = theIonTable->GetParticle(i)->GetProcessManager();
	  pmanager->SetVerboseLevel(3);
	  pmanager ->AddProcess(theRadioactiveDecay);
	  pmanager ->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
	  pmanager ->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
	}
    } 
}

void Tst18PhysicsList::SetCuts()
{
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}

void Tst18PhysicsList::ConstructAllBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst18PhysicsList::ConstructAllLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst18PhysicsList::ConstructAllMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst18PhysicsList::ConstructAllBaryons()
{
  //  Construct all barions
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst18PhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst18PhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

