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

#include <iomanip>                

#include "globals.hh"
#include "Tst13PhysicsList.hh"
#include "G4SystemOfUnits.hh"
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

Tst13PhysicsList::Tst13PhysicsList():  G4VUserPhysicsList()
{
  SetVerboseLevel(1);
}

Tst13PhysicsList::~Tst13PhysicsList()
{}

void Tst13PhysicsList::ConstructParticle()
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

void Tst13PhysicsList::ConstructAllBosons()
{
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst13PhysicsList::ConstructAllLeptons()
{
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst13PhysicsList::ConstructAllMesons()
{
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst13PhysicsList::ConstructAllBaryons()
{
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst13PhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst13PhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst13PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructLeptHad();
  ConstructHad();
  ConstructGeneral();
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

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

void Tst13PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
      // Construct processes for gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
      // Construct processes for electron
      pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4eIonisation(),-1,2,2);
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);
  
    } else if (particleName == "e+") {
      // Construct processes for positron
      pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4eIonisation(),-1,2,2);
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);      
      pmanager->AddProcess(new G4eplusAnnihilation(),0,-1,4);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
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
          (particle->GetParticleName() != "chargedgeantino") &&
          (!particle->IsShortLived()) ) {  
     // all others charged particles except geantino
       pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
       pmanager->AddProcess(new G4hIonisation(),-1,2,2);       
      }
    }
  }
}

// Hadron Processes

#include "G4HadronElasticProcess.hh"
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

// Elastic scattering model
#include "G4HadronElastic.hh"

// Medium-energy models
#include "G4CascadeInterface.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4QMDReaction.hh"
#include "G4NeutronRadCapture.hh"

// Neutron capture cross sections
#include "G4NeutronCaptureXS.hh"

// High-energy and generator models
#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4Evaporation.hh"
#include "G4CompetitiveFission.hh"
#include "G4FermiBreakUp.hh"
#include "G4StatMF.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4StringModel.hh"
#include "G4PreCompoundModel.hh"
#include "G4QGSModel.hh"
#include "G4FTFModel.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSMFragmentation.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
 
// ConstructHad()
//
// Makes discrete physics processes for the hadrons
// The processes are: Elastic scattering, Inelastic scattering,
// Fission (for neutron only), and Capture (neutron).
//

void Tst13PhysicsList::ConstructHad()
{
  // Models for low and medium energy hadrons
  G4CascadeInterface* bertini = new G4CascadeInterface;
  bertini->SetMinEnergy(0.0);
  bertini->SetMaxEnergy(20.0*GeV);

  // Models for low and medium energy ions
  G4BinaryLightIonReaction* binaryLightIon = new G4BinaryLightIonReaction;
  binaryLightIon->SetMinEnergy(0.0);
  binaryLightIon->SetMaxEnergy(110*MeV);

  G4QMDReaction* qmd = new G4QMDReaction;
  qmd->SetMinEnergy(100*MeV);
  qmd->SetMaxEnergy(10*GeV);

  // Model for high energy nucleons: QGSP for protons and neutrons above 19 GeV
  G4TheoFSGenerator* theTheoModel = new G4TheoFSGenerator;
       
  // all models for treatment of thermal nucleus 
  G4Evaporation* theEvaporation = new G4Evaporation;
  G4FermiBreakUp* theFermiBreakUp = new G4FermiBreakUp;
  G4StatMF* theMF = new G4StatMF;

  // Evaporation logic
  G4ExcitationHandler* theHandler = new G4ExcitationHandler;
  theHandler->SetEvaporation(theEvaporation);
  theHandler->SetFermiModel(theFermiBreakUp);
  theHandler->SetMultiFragmentation(theMF);
  theHandler->SetMaxAandZForFermiBreakUp(12, 6);
  theHandler->SetMinEForMultiFrag(10*MeV);
	
  // Pre equilibrium stage (dummy for now)
  G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel(theHandler);
    
  // a no-cascade generator-precompound interface
  G4GeneratorPrecompoundInterface* theCascade = new G4GeneratorPrecompoundInterface;
  theCascade->SetDeExcitation(thePreEquilib);  
	
  // here come the high energy parts
  G4VPartonStringModel* theStringModel;
  theStringModel = new G4QGSModel<G4QGSParticipants>;
  theTheoModel->SetTransport(theCascade);
  theTheoModel->SetHighEnergyGenerator(theStringModel);
  theTheoModel->SetMinEnergy(19*GeV);
  theTheoModel->SetMaxEnergy(100*TeV);

  G4VLongitudinalStringDecay* theFragmentation = new G4QGSMFragmentation;
  G4ExcitedStringDecay* theStringDecay = new G4ExcitedStringDecay(theFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);
  // done with the generator model (most of the above is also available as default)

  // Basic elastic scattering - default energy range and cross sections
  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  G4HadronElastic* elasticModel = new G4HadronElastic();
  theElasticProcess->RegisterMe(elasticModel);

  // FTFP model for hadrons above 19 GeV 
  G4TheoFSGenerator* ftfp_hi = new G4TheoFSGenerator;
  G4FTFModel* stringModel = new G4FTFModel;
  G4ExcitedStringDecay* stringDecay =
    new G4ExcitedStringDecay(new G4LundStringFragmentation);
  stringModel->SetFragmentationModel(stringDecay);

  G4GeneratorPrecompoundInterface* cascade = new G4GeneratorPrecompoundInterface;
  G4PreCompoundModel* preEquilib = new G4PreCompoundModel(new G4ExcitationHandler);
  cascade->SetDeExcitation(preEquilib);  

  ftfp_hi->SetTransport(cascade);
  ftfp_hi->SetHighEnergyGenerator(stringModel);
  ftfp_hi->SetMinEnergy(19*GeV);
  ftfp_hi->SetMaxEnergy(100*TeV);

  // FTFP model for anti-baryons of all energies
  G4TheoFSGenerator* ftfp_anti = new G4TheoFSGenerator;
  ftfp_anti->SetTransport(cascade);
  ftfp_anti->SetHighEnergyGenerator(stringModel);
  ftfp_anti->SetMinEnergy(0.0);
  ftfp_anti->SetMaxEnergy(100*TeV);

  theParticleIterator->reset();
  while ((*theParticleIterator)()) {
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "pi+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4PionPlusInelasticProcess* theInelasticProcess = 
                             new G4PionPlusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(ftfp_hi);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "pi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4PionMinusInelasticProcess* theInelasticProcess = 
                             new G4PionMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(ftfp_hi);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonPlusInelasticProcess* theInelasticProcess = 
                             new G4KaonPlusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(ftfp_hi);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon0S") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonZeroSInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroSInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(ftfp_hi);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon0L") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonZeroLInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroLInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(ftfp_hi);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "kaon-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonMinusInelasticProcess* theInelasticProcess = 
                             new G4KaonMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(ftfp_hi);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "proton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4ProtonInelasticProcess* theInelasticProcess = 
                             new G4ProtonInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_proton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiProtonInelasticProcess* theInelasticProcess = 
                             new G4AntiProtonInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(ftfp_anti);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "neutron") {     
      // elastic scattering
      pmanager->AddDiscreteProcess(theElasticProcess);

      // inelastic scattering
      G4NeutronInelasticProcess* theInelasticProcess = 
                             new G4NeutronInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);

      // capture
      G4HadronCaptureProcess* theCaptureProcess = new G4HadronCaptureProcess;
      G4NeutronRadCapture* captureModel = new G4NeutronRadCapture;
      theCaptureProcess->RegisterMe(captureModel);
      G4NeutronCaptureXS* theCaptureXS = new G4NeutronCaptureXS;
      theCaptureProcess->AddDataSet(theCaptureXS);
      pmanager->AddDiscreteProcess(theCaptureProcess);

    } else if (particleName == "anti_neutron") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiNeutronInelasticProcess* theInelasticProcess = 
                             new G4AntiNeutronInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(ftfp_anti);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "lambda") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LambdaInelasticProcess* theInelasticProcess = 
                             new G4LambdaInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(ftfp_hi);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_lambda") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiLambdaInelasticProcess* theInelasticProcess = 
                             new G4AntiLambdaInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(ftfp_anti);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "sigma+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4SigmaPlusInelasticProcess* theInelasticProcess = 
                             new G4SigmaPlusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(ftfp_hi);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "sigma-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4SigmaMinusInelasticProcess* theInelasticProcess = 
                             new G4SigmaMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(ftfp_hi);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_sigma+") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiSigmaPlusInelasticProcess* theInelasticProcess = 
                             new G4AntiSigmaPlusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(ftfp_anti);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_sigma-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiSigmaMinusInelasticProcess* theInelasticProcess = 
                             new G4AntiSigmaMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(ftfp_anti);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "xi0") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4XiZeroInelasticProcess* theInelasticProcess = 
                             new G4XiZeroInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(ftfp_hi);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "xi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4XiMinusInelasticProcess* theInelasticProcess = 
                             new G4XiMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(ftfp_hi);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "omega-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4OmegaMinusInelasticProcess* theInelasticProcess =
                             new G4OmegaMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(bertini);
      theInelasticProcess->RegisterMe(ftfp_hi);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_xi0") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiXiZeroInelasticProcess* theInelasticProcess = 
                             new G4AntiXiZeroInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(ftfp_anti);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_xi-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiXiMinusInelasticProcess* theInelasticProcess = 
                             new G4AntiXiMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(ftfp_anti);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "anti_omega-") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiOmegaMinusInelasticProcess* theInelasticProcess =
                             new G4AntiOmegaMinusInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(ftfp_anti);
      pmanager->AddDiscreteProcess(theInelasticProcess);
 
    } else if (particleName == "deuteron") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4DeuteronInelasticProcess* theInelasticProcess = 
                             new G4DeuteronInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(binaryLightIon);
      theInelasticProcess->RegisterMe(qmd);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "triton") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4TritonInelasticProcess* theInelasticProcess = 
                             new G4TritonInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(binaryLightIon);
      theInelasticProcess->RegisterMe(qmd);
      pmanager->AddDiscreteProcess(theInelasticProcess);

    } else if (particleName == "alpha") {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AlphaInelasticProcess* theInelasticProcess = 
                             new G4AlphaInelasticProcess("inelastic");
      theInelasticProcess->RegisterMe(binaryLightIon);
      theInelasticProcess->RegisterMe(qmd);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  } // while
}

void Tst13PhysicsList::ConstructLeptHad()
{}

#include "G4Decay.hh"
void Tst13PhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}


void Tst13PhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "Tst13PhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue/mm << " (mm)" << G4endl;
  }  
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}


