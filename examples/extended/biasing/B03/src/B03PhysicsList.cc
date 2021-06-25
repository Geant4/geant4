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
/// \file biasing/B03/src/B03PhysicsList.cc
/// \brief Implementation of the B03PhysicsList class
//
//
//

#include "globals.hh"
#include <iomanip>                

#include "B03PhysicsList.hh"

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
#include "G4SystemOfUnits.hh"
#include "G4HadronicParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B03PhysicsList::B03PhysicsList(G4String parallelname):  
G4VUserPhysicsList(),
fBiasWorldName(parallelname)
{
  fParaWorldName.clear();
  SetVerboseLevel(1);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B03PhysicsList::~B03PhysicsList()
{
  fParaWorldName.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B03PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructAllBosons();
  ConstructAllLeptons();
  ConstructAllMesons();
  ConstructAllBaryons();
  ConstructAllIons();
  ConstructAllShortLiveds();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B03PhysicsList::ConstructAllBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B03PhysicsList::ConstructAllLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B03PhysicsList::ConstructAllMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B03PhysicsList::ConstructAllBaryons()
{
  //  Construct all barions
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B03PhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B03PhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B03PhysicsList::ConstructProcess()
{
  AddTransportation();
  AddScoringProcess();
  AddBiasingProcess();
  ConstructEM();
  ConstructLeptHad();
  ConstructHad();
  ConstructGeneral();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

void B03PhysicsList::ConstructEM()
{
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
      // Construct processes for gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
    //electron
      // Construct processes for electron
      pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4eIonisation(),-1,2,2);
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);
  
    } else if (particleName == "e+") {
    //positron
      // Construct processes for positron
     pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
     
     pmanager->AddProcess(new G4eIonisation(),-1,2,2);
     pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);      
     pmanager->AddProcess(new G4eplusAnnihilation(),0,-1,4);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
    //muon  
     // Construct processes for muon+
     pmanager->AddProcess(new G4MuMultipleScattering(),-1,1,1);
     pmanager->AddProcess(new G4MuIonisation(),-1,2,2);
     pmanager->AddProcess(new G4MuBremsstrahlung(),-1,-1,3);
     pmanager->AddProcess(new G4MuPairProduction(),-1,-1,4);       
     
    } else if( particleName == "GenericIon" ) {
      pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4hIonisation(),-1,2,2); 
    } else { 
      if ((particle->GetPDGCharge() != 0.0) && 
          (particle->GetParticleName() != "chargedgeantino")&&
          (!particle->IsShortLived()) ) {
     // all others charged particles except geantino
       pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
       pmanager->AddProcess(new G4hIonisation(),-1,2,2);       
     }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Hadron Processes

#include "G4HadronElasticProcess.hh"
#include "G4NeutronFissionProcess.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4HadronInelasticProcess.hh"

// Low-energy Models

#include "G4HadronElastic.hh"
#include "G4NeutronRadCapture.hh"
#include "G4LFission.hh"

// -- generator models
#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4CompetitiveFission.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4CascadeInterface.hh"
#include "G4StringModel.hh"
#include "G4PreCompoundModel.hh"
#include "G4FTFModel.hh"
#include "G4QGSMFragmentation.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QMDReaction.hh"
#include "G4BinaryLightIonReaction.hh"

// Cross sections
#include "G4ComponentGGNuclNuclXsc.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4CrossSectionElastic.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4NeutronInelasticXS.hh"

//
// ConstructHad()
//
// Makes discrete physics processes for the hadrons
// The processes are: Elastic scattering, Inelastic scattering,
// Fission (for neutron only), and Capture (neutron).
//

void B03PhysicsList::ConstructHad()
{
  // this will be the model class for high energies
  G4TheoFSGenerator* theTheoModel = new G4TheoFSGenerator;
  G4TheoFSGenerator* antiBHighEnergyModel = new G4TheoFSGenerator;
       
  // Evaporation logic
  G4ExcitationHandler* theHandler = new G4ExcitationHandler;
  theHandler->SetMinEForMultiFrag(3*MeV);
        
  // Pre equilibrium stage 
  G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel(theHandler);

  // a no-cascade generator-precompound interaface
  G4GeneratorPrecompoundInterface* theCascade =
                                    new G4GeneratorPrecompoundInterface;
  theCascade->SetDeExcitation(thePreEquilib);  

  // Bertini cascade
  G4CascadeInterface* bertini = new G4CascadeInterface;
  bertini->SetMaxEnergy(22*MeV);

  // here come the high energy parts
  G4VPartonStringModel* theStringModel;
  theStringModel = new G4FTFModel;
  theTheoModel->SetTransport(theCascade);
  theTheoModel->SetHighEnergyGenerator(theStringModel);
  theTheoModel->SetMinEnergy(19*GeV);
  theTheoModel->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );

  G4VLongitudinalStringDecay* theFragmentation = new G4QGSMFragmentation;
  G4ExcitedStringDecay* theStringDecay = 
                         new G4ExcitedStringDecay(theFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);

  // high energy model for anti-baryons
  antiBHighEnergyModel = new G4TheoFSGenerator("ANTI-FTFP");
  G4FTFModel* antiBStringModel = new G4FTFModel;
  G4ExcitedStringDecay* stringDecay = 
                    new G4ExcitedStringDecay(new G4LundStringFragmentation);
  antiBStringModel->SetFragmentationModel(stringDecay);

  G4GeneratorPrecompoundInterface* antiBCascade = 
                               new G4GeneratorPrecompoundInterface;
  G4PreCompoundModel* preEquilib = 
                  new G4PreCompoundModel(new G4ExcitationHandler);
  antiBCascade->SetDeExcitation(preEquilib);

  antiBHighEnergyModel->SetTransport(antiBCascade);
  antiBHighEnergyModel->SetHighEnergyGenerator(antiBStringModel);
  antiBHighEnergyModel->SetMinEnergy(0.0);
  antiBHighEnergyModel->SetMaxEnergy(20*TeV);

  // Light ion models
  G4BinaryLightIonReaction* binaryCascade = new G4BinaryLightIonReaction;
  binaryCascade->SetMinEnergy(0.0);
  binaryCascade->SetMaxEnergy(110*MeV);

  G4QMDReaction* qmd = new G4QMDReaction;
  qmd->SetMinEnergy(100*MeV);
  qmd->SetMaxEnergy(10*GeV);

  G4VCrossSectionDataSet* ionXS = new G4CrossSectionInelastic( new G4ComponentGGNuclNuclXsc );

  G4ComponentGGHadronNucleusXsc * ggHNXsec = new G4ComponentGGHadronNucleusXsc();
  G4VCrossSectionDataSet * theGGHNEl = new G4CrossSectionElastic(ggHNXsec);
  G4VCrossSectionDataSet * theGGHNInel = new G4CrossSectionInelastic(ggHNXsec);

  // Elastic process
  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  theElasticProcess->AddDataSet(theGGHNEl);
  G4HadronElastic* theElasticModel = new G4HadronElastic;
  theElasticProcess->RegisterMe(theElasticModel);

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)()) {
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "pi+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4PionPlus::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(bertini);
         theInelasticProcess->RegisterMe(theTheoModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
    } else if (particleName == "pi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4PionMinus::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(bertini);
         theInelasticProcess->RegisterMe(theTheoModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
    } else if (particleName == "kaon+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4KaonPlus::Definition() );	 
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(bertini);
         theInelasticProcess->RegisterMe(theTheoModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon0S") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4KaonZeroShort::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(bertini);
         theInelasticProcess->RegisterMe(theTheoModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon0L") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4KaonZeroLong::Definition() );	 
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(bertini);
         theInelasticProcess->RegisterMe(theTheoModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4KaonMinus::Definition() );	 
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(bertini);
         theInelasticProcess->RegisterMe(theTheoModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Proton::Definition() );	 
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(bertini);
         theInelasticProcess->RegisterMe(theTheoModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiProton::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(antiBHighEnergyModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "neutron") {         
      // elastic scattering
      pmanager->AddDiscreteProcess(theElasticProcess);

      // inelastic scattering
      G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Neutron::Definition() );
      theInelasticProcess->AddDataSet(new G4NeutronInelasticXS());
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

      // fission
      G4NeutronFissionProcess* theFissionProcess = new G4NeutronFissionProcess;
      G4LFission* theFissionModel = new G4LFission;
      theFissionProcess->RegisterMe(theFissionModel);
      pmanager->AddDiscreteProcess(theFissionProcess);

      // capture
      G4NeutronCaptureProcess* theCaptureProcess = new G4NeutronCaptureProcess;
      G4NeutronRadCapture* theCaptureModel = new G4NeutronRadCapture;
      theCaptureProcess->RegisterMe(theCaptureModel);
      pmanager->AddDiscreteProcess(theCaptureProcess);

    } else if (particleName == "anti_neutron") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiNeutron::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(antiBHighEnergyModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "lambda") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Lambda::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(bertini);
         theInelasticProcess->RegisterMe(theTheoModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_lambda") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiLambda::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(antiBHighEnergyModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "sigma+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4SigmaPlus::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(bertini);
         theInelasticProcess->RegisterMe(theTheoModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "sigma-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4SigmaMinus::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(bertini);
         theInelasticProcess->RegisterMe(theTheoModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_sigma+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiSigmaPlus::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(antiBHighEnergyModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_sigma-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiSigmaMinus::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(antiBHighEnergyModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "xi0") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4XiZero::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(bertini);
         theInelasticProcess->RegisterMe(theTheoModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "xi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4XiMinus::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(bertini);
         theInelasticProcess->RegisterMe(theTheoModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_xi0") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiXiZero::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(antiBHighEnergyModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_xi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4AntiXiMinus::Definition() );
         theInelasticProcess->AddDataSet(theGGHNInel);
         theInelasticProcess->RegisterMe(antiBHighEnergyModel);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "deuteron") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Deuteron::Definition() );
         theInelasticProcess->RegisterMe(binaryCascade);
         theInelasticProcess->RegisterMe(qmd);
         theInelasticProcess->AddDataSet(ionXS);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "triton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Triton::Definition() );
         theInelasticProcess->RegisterMe(binaryCascade);
         theInelasticProcess->RegisterMe(qmd);
         theInelasticProcess->AddDataSet(ionXS);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "alpha") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4HadronInelasticProcess* theInelasticProcess = 
	   new G4HadronInelasticProcess( "inelastic", G4Alpha::Definition() );
         theInelasticProcess->RegisterMe(binaryCascade);
         theInelasticProcess->RegisterMe(qmd);
         theInelasticProcess->AddDataSet(ionXS);
         pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "omega-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
        new G4HadronInelasticProcess( "inelastic", G4OmegaMinus::Definition() );
      theInelasticProcess->AddDataSet(theGGHNInel);
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_omega-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4HadronInelasticProcess* theInelasticProcess = 
        new G4HadronInelasticProcess( "inelastic", G4AntiOmegaMinus::Definition() );
      theInelasticProcess->AddDataSet(theGGHNInel);
      theInelasticProcess->RegisterMe(antiBHighEnergyModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B03PhysicsList::ConstructLeptHad()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Decay.hh"
void B03PhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B03PhysicsList::SetCuts()
{
  if (verboseLevel >0)
  {
    G4cout << "B03PhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue/mm << " (mm)" << G4endl;
  }  
  //   "G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ParallelWorldProcess.hh"
void B03PhysicsList::AddScoringProcess(){

  G4int npw = fParaWorldName.size();
  for ( G4int i = 0; i < npw; i++){
   G4String procName = "ParaWorldProc_"+fParaWorldName[i];
   G4ParallelWorldProcess* theParallelWorldProcess
     = new G4ParallelWorldProcess(procName);
   theParallelWorldProcess->SetParallelWorld(fParaWorldName[i]);

   auto particleIterator=GetParticleIterator();
   particleIterator->reset();
   while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    pmanager->AddProcess(theParallelWorldProcess);
    if(theParallelWorldProcess->IsAtRestRequired(particle))
    {pmanager->SetProcessOrdering(theParallelWorldProcess, idxAtRest, 9900);}
    pmanager->SetProcessOrderingToSecond(theParallelWorldProcess, idxAlongStep);
    pmanager->SetProcessOrdering(theParallelWorldProcess, idxPostStep, 9900);
   }
  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ImportanceProcess.hh"
#include "G4IStore.hh"
void B03PhysicsList::AddBiasingProcess(){

  G4cout << " Preparing Importance Sampling with GhostWorld " 
         << fBiasWorldName << G4endl;
  
  G4IStore* iStore = G4IStore::GetInstance(fBiasWorldName);
  G4GeometrySampler fGeomSampler(fBiasWorldName,"neutron");
  fGeomSampler.SetParallel(true); // parallelworld
  // fGeomSampler.SetWorld(iStore->GetParallelWorldVolumePointer());
  //  fGeomSampler->PrepareImportanceSampling(G4IStore::
  //                              GetInstance(fBiasWorldName), 0);
  static G4bool first = true;
  if(first) {
    fGeomSampler.PrepareImportanceSampling(iStore, 0);

    fGeomSampler.Configure();
    G4cout << " GeomSampler Configured!!! " << G4endl;
    first = false;
  }

#ifdef G4MULTITHREADED 
  if(!G4Threading::IsMasterThread()) fGeomSampler.AddProcess();
#else
  G4cout << " Running in singlethreaded mode!!! " << G4endl;
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
