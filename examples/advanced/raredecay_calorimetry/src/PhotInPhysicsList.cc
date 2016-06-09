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
// $Id: PhotInPhysicsList.cc,v 1.4 2006/06/29 16:25:19 gunter Exp $
// GEANT4 tag $Name: geant4-08-02 $
//

#define debug

#include "PhotInPhysicsList.hh"

PhotInPhysicsList::PhotInPhysicsList():  G4VUserPhysicsList() { SetVerboseLevel(1); }

PhotInPhysicsList::~PhotInPhysicsList() {}

void PhotInPhysicsList::ConstructParticle()
{
  // In this method, static member functions for particles should be called
  // for all particles which user is going to use in the simulation. If not
  // defined particle appear in the simulation, it can cause a WORNING which
  // means that in the simulation appeard unexpected particles. Then add them.

  //// @@ Word "Definition" can be skipped. - Old fashion(M.K.)
  //
  //// pseudo-particles
  //G4Geantino::GeantinoDefinition();
  //G4ChargedGeantino::ChargedGeantinoDefinition();
  //
  //// gammas
  //G4Gamma::GammaDefinition();
  //
  //// leptons (without tau and it's neutrino)
  //G4Electron::ElectronDefinition();
  //G4Positron::PositronDefinition();
  //G4MuonPlus::MuonPlusDefinition();
  //G4MuonMinus::MuonMinusDefinition();
  //
  //G4NeutrinoE::NeutrinoEDefinition();
  //G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  //G4NeutrinoMu::NeutrinoMuDefinition();
  //G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
  //
  ////  mesons
  //G4PionPlus::PionPlusDefinition();
  //G4PionMinus::PionMinusDefinition();
  //G4PionZero::PionZeroDefinition();
  //G4Eta::EtaDefinition();
  //G4EtaPrime::EtaPrimeDefinition();
  //G4KaonPlus::KaonPlusDefinition();
  //G4KaonMinus::KaonMinusDefinition();
  //G4KaonZero::KaonZeroDefinition();
  //G4AntiKaonZero::AntiKaonZeroDefinition();
  //G4KaonZeroLong::KaonZeroLongDefinition();
  //G4KaonZeroShort::KaonZeroShortDefinition();
  //
  ////  barions
  //G4Proton::ProtonDefinition();
  //G4AntiProton::AntiProtonDefinition();
  //G4Neutron::NeutronDefinition();
  //G4AntiNeutron::AntiNeutronDefinition();
  //
  //// hyperons
		//G4Lambda::LambdaDefinition();
		//G4SigmaPlus::SigmaPlusDefinition();
		//G4SigmaZero::SigmaZeroDefinition();
		//G4SigmaMinus::SigmaMinusDefinition();
		//G4XiMinus::XiMinusDefinition();
		//G4XiZero::XiZeroDefinition();
		//G4OmegaMinus::OmegaMinusDefinition();
		//G4AntiLambda::AntiLambdaDefinition();
		//G4AntiSigmaPlus::AntiSigmaPlusDefinition();
		//G4AntiSigmaZero::AntiSigmaZeroDefinition();
		//G4AntiSigmaMinus::AntiSigmaMinusDefinition();
		//G4AntiXiMinus::AntiXiMinusDefinition();
		//G4AntiXiZero::AntiXiZeroDefinition();
		//G4AntiOmegaMinus::AntiOmegaMinusDefinition();

  // The same can be done shorter, using supported Geant4 particle initialization classes:
		G4BosonConstructor      consBos; consBos.ConstructParticle();
  G4LeptonConstructor     consLep; consLep.ConstructParticle();
  G4MesonConstructor      consMes; consMes.ConstructParticle();
  G4BaryonConstructor     consBar; consBar.ConstructParticle();
  G4IonConstructor        consIon; consIon.ConstructParticle();
  G4ShortLivedConstructor consShL; consShL.ConstructParticle();
}


void PhotInPhysicsList::ConstructProcess()
{
#ifdef debug
  G4cout<<"PhotInPhysicsList::ConstructProcess: is called "<<G4endl;
#endif
  AddTransportation();     // Transportation is a "process" and defined in the basic class

  // Add Electromagnetic interaction Processes and Decays
  G4Decay* theDecayProcess = new G4Decay(); // @@ When this class is decayed? (M.K.)

  // For high energy hadronic (QGSC): this will be the model class for high energies
  G4TheoFSGenerator* theTheoModel = new G4TheoFSGenerator;
         
  // a cascade interface interface to chiral invriant phase space decay
  // at particle level.
  G4StringChipsParticleLevelInterface* theCascade =
                                                 new G4StringChipsParticleLevelInterface();
	
  // here come the high energy parts
  // the string model; still not quite according to design - Explicite use of the forseen
  // interfaces  will be tested and documented in this program by beta-02 at latest.
  G4VPartonStringModel* theStringModel;
  theStringModel = new G4QGSModel<G4QGSParticipants>;
  theTheoModel->SetTransport(theCascade);
  theTheoModel->SetHighEnergyGenerator(theStringModel);
  theTheoModel->SetMinEnergy(19*GeV);
  theTheoModel->SetMaxEnergy(100*TeV);

  G4VLongitudinalStringDecay* theFragmentation = new G4LundStringFragmentation;
  G4ExcitedStringDecay* theStringDecay = new G4ExcitedStringDecay(theFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);

  // done with the generator model (most of the above is also available as default)
  G4HadronElasticProcess* theElasticProcess =  new G4HadronElasticProcess; //@@
  G4LElastic* theElasticModel = new G4LElastic;                            //@@
  theElasticProcess->RegisterMe(theElasticModel);
  G4HadronElasticProcess* theElasticProcess1 =  new G4HadronElasticProcess;

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
  {
    G4ParticleDefinition* particle = theParticleIterator->value();
#ifdef debug
				G4cout<<"PhotInPhysList::ConstructProcess: Part="<<particle->GetParticleName()<<G4endl;
#endif
    G4ProcessManager* pmanager = particle->GetProcessManager();

    // ================================ Decays ===========================================

    if (theDecayProcess->IsApplicable(*particle))
    { 
      pmanager ->AddProcess(theDecayProcess);
      // set ordering for PostStepDoIt and AtRestDoIt
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }

    G4String particleName = particle->GetParticleName();

    // =========================== EM Interactions & Photo-nuclear =======================

    if (particleName == "gamma")    // gamma
    {
      // Construct processes for gamma
      // ........The OLD (GHAR) usage of the CHIPS Photo-nuclear reactions .............
      //G4PhotoNuclearProcess* thePhotoNuclearProcess = new G4PhotoNuclearProcess;
      //G4QGSMFragmentation* theFragmentation = new G4QGSMFragmentation;
						//G4QGSModel<G4GammaParticipants>*theStringModel=new G4QGSModel<G4GammaParticipants>;
      //G4GammaNuclearReaction* theGammaReaction = new G4GammaNuclearReaction;
      //G4TheoFSGenerator* theModel = new G4TheoFSGenerator;
      //G4StringChipsParticleLevelInterface* theCascade =
      //                                           new G4StringChipsParticleLevelInterface;
						//theModel->SetTransport(theCascade);
      //theModel->SetHighEnergyGenerator(theStringModel);
      //G4ExcitedStringDecay* theStringDecay = new G4ExcitedStringDecay(theFragmentation);
      //theStringModel->SetFragmentationModel(theStringDecay);
      //theGammaReaction->SetMaxEnergy(3.5*GeV);
      //thePhotoNuclearProcess->RegisterMe(theGammaReaction);
      //theModel->SetMinEnergy(3.*GeV);
      //theModel->SetMaxEnergy(100*TeV);
      //thePhotoNuclearProcess->RegisterMe(theModel);
      //pmanager->AddDiscreteProcess(thePhotoNuclearProcess);
      // ........The NEW (native) usage of the CHIPS Photo-nuclear reactions .............
      pmanager->AddDiscreteProcess(new G4QCollision);

      // Pure Electromagnetic Processes
      pmanager->AddDiscreteProcess(new G4GammaConversion());     // Pair Production
      pmanager->AddDiscreteProcess(new G4ComptonScattering());   // Compton Effect     
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect()); // Photo Effect
    } 
    else if (particleName == "e-")    //electron
    {
      // Construct processes for electron
      // ........The OLD (GHAR) usage of the CHIPS Photo-nuclear reactions .............
      //G4ElectronNuclearProcess* theElectronNuclearProcess = new G4ElectronNuclearProcess;
      //G4ElectroNuclearReaction* theElectronReaction = new G4ElectroNuclearReaction;
      //theElectronNuclearProcess->RegisterMe(theElectronReaction);
      ////theElectronNuclearProcess->BiasCrossSectionByFactor(100); // Very dangerous
      //pmanager->AddDiscreteProcess(theElectronNuclearProcess);
      // ........The NEW (native) usage of the CHIPS electro-nuclear reactions ............
      pmanager->AddDiscreteProcess(new G4QCollision);

      // Pure Electromagnetic Processes
      pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);// ElectronMultipleScattering
      pmanager->AddProcess(new G4eIonisation(),-1,2,2);       // Electron Ionisation							
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);  // Electron BremsStrahlung   
    }
    else if (particleName == "e+")    //positron
    {
      // Construct processes for positron
      // ........The OLD (GHAR) usage of the CHIPS Photo-nuclear reactions .............
      //G4PositronNuclearProcess* thePositronNuclearProcess = new G4PositronNuclearProcess;
      //G4ElectroNuclearReaction* thePositronReaction = new G4ElectroNuclearReaction;
      //thePositronNuclearProcess->RegisterMe(thePositronReaction);
      ////thePositronNuclearProcess->BiasCrossSectionByFactor(100); // Very dangerous
      //pmanager->AddDiscreteProcess(thePositronNuclearProcess);
      // ........The NEW (native) usage of the CHIPS positron-nuclear reactions ...........
      pmanager->AddDiscreteProcess(new G4QCollision);

      // Pure Electromagnetic Processes
      pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);// PositronMultipleScattering
      pmanager->AddProcess(new G4eIonisation(),-1,2,2);							// Positron Ionisation
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);  // Positron BremsStrahlung
      G4eplusAnnihilation* theAnnihilation = new G4eplusAnnihilation;
      pmanager->AddDiscreteProcess(theAnnihilation);     // Positron Annihilation on Flight
      pmanager->AddRestProcess(theAnnihilation);         // Positron Annihilation at Rest  
    }
    else if( particleName == "mu+" || particleName == "mu-" )  //muon  of both signs
    {
      // Construct processes for muon+/-
      pmanager->AddProcess(new G4MultipleScattering(),-1,1,1); // Mu Multiple Scattering
      pmanager->AddProcess(new G4MuIonisation(),-1,2,2);       // Mu Ionization
      pmanager->AddProcess(new G4MuBremsstrahlung(),-1,-1,3);  // Mu Bremsstrahlung
      pmanager->AddProcess(new G4MuPairProduction(),-1,-1,4);  // Mu Pair Production      
      // ........The NEW (native) usage of the CHIPS Muo-nuclear reactions .............
      pmanager->AddDiscreteProcess(new G4QCollision);
      // ........Negative lepton stopping
      if(particleName == "mu-") pmanager->AddRestProcess(new G4QCaptureAtRest, ordDefault);
    }
    else if( particleName == "tau+" || particleName == "tau-" )
    {
      // Construct processes for muon+/-
      pmanager->AddProcess(new G4MultipleScattering(),-1,1,1); // Tau Multiple Scattering
      // ........The NEW (native) usage of the CHIPS Muo-nuclear reactions .............
      pmanager->AddDiscreteProcess(new G4QCollision);
      // ........Negative lepton stopping
      if(particleName == "tau-")pmanager->AddRestProcess(new G4QCaptureAtRest, ordDefault);
    }
    else if( particleName == "GenericIon" )
    {
      // Construct processes for ions
      pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4hIonisation(),-1,2,2); 
    }
    else
    { 
      if (particle->GetPDGCharge() && (particle->GetParticleName() != "chargedgeantino") &&
                                                                !particle->IsShortLived() )
      {  
        // short lived particles except geantino
        pmanager->AddProcess(new G4MultipleScattering(),-1,1,1);
        pmanager->AddProcess(new G4hIonisation(),-1,2,2);       
      }
    }

    // ============================== Hadronic Interactions ===============================

    if (particleName == "pi+")
    {
       pmanager->AddDiscreteProcess(theElasticProcess);
       G4PionPlusInelasticProcess* theInelasticProcess = 
                                new G4PionPlusInelasticProcess("inelastic");
       G4LEPionPlusInelastic* theInelasticModel =  new G4LEPionPlusInelastic;
       theInelasticProcess->RegisterMe(theInelasticModel);
       theInelasticProcess->RegisterMe(theTheoModel);
       pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "pi-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4PionMinusInelasticProcess* theInelasticProcess = 
                                new G4PionMinusInelasticProcess("inelastic");
      G4LEPionMinusInelastic* theInelasticModel =  new G4LEPionMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4QCaptureAtRest, ordDefault);
    }
    else if (particleName == "kaon+")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonPlusInelasticProcess* theInelasticProcess = 
                                  new G4KaonPlusInelasticProcess("inelastic");
      G4LEKaonPlusInelastic* theInelasticModel = new G4LEKaonPlusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      new G4HEKaonPlusInelastic; // @@
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "kaon0S")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonZeroSInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroSInelasticProcess("inelastic");
      G4LEKaonZeroSInelastic* theInelasticModel =  new G4LEKaonZeroSInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "kaon0L")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonZeroLInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroLInelasticProcess("inelastic");
      G4LEKaonZeroLInelastic* theInelasticModel =  new G4LEKaonZeroLInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      new G4HEKaonZeroInelastic;
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "kaon-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4KaonMinusInelasticProcess* theInelasticProcess = 
                                 new G4KaonMinusInelasticProcess("inelastic");
      G4LEKaonMinusInelastic* theInelasticModel =  new G4LEKaonMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4QCaptureAtRest, ordDefault);
    }
    else if (particleName == "proton")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4ProtonInelasticProcess* theInelasticProcess = 
                                    new G4ProtonInelasticProcess("inelastic");
      G4LEProtonInelastic* theInelasticModel = new G4LEProtonInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_proton")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiProtonInelasticProcess* theInelasticProcess = 
                                new G4AntiProtonInelasticProcess("inelastic");
      G4LEAntiProtonInelastic* theInelasticModel =  new G4LEAntiProtonInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiProtonInelastic* theHEInelasticModel =  new G4HEAntiProtonInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4QCaptureAtRest, ordDefault);
    }
    else if (particleName == "neutron")
    {     
      // elastic scattering
      G4LElastic* theElasticModel1 = new G4LElastic;
      theElasticProcess1->RegisterMe(theElasticModel1);
      pmanager->AddDiscreteProcess(theElasticProcess1);
      // inelastic scattering
      G4NeutronInelasticProcess* theInelasticProcess = 
                                    new G4NeutronInelasticProcess("inelastic");
      G4LENeutronInelastic* theInelasticModel = new G4LENeutronInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      theInelasticProcess->RegisterMe(theTheoModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      // fission
      G4HadronFissionProcess* theFissionProcess = new G4HadronFissionProcess;
      G4LFission* theFissionModel = new G4LFission;
      theFissionProcess->RegisterMe(theFissionModel);
      pmanager->AddDiscreteProcess(theFissionProcess);
      // capture
      G4HadronCaptureProcess* theCaptureProcess = new G4HadronCaptureProcess;
      G4LCapture* theCaptureModel = new G4LCapture;
      theCaptureProcess->RegisterMe(theCaptureModel);
      pmanager->AddDiscreteProcess(theCaptureProcess);
    }  
    else if (particleName == "anti_neutron")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiNeutronInelasticProcess* theInelasticProcess = 
                               new G4AntiNeutronInelasticProcess("inelastic");
      G4LEAntiNeutronInelastic* theInelasticModel =  new G4LEAntiNeutronInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiNeutronInelastic* theHEInelasticModel =  new G4HEAntiNeutronInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4QCaptureAtRest, ordDefault);
    }
    else if (particleName == "lambda")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4LambdaInelasticProcess* theInelasticProcess = 
                                    new G4LambdaInelasticProcess("inelastic");
      G4LELambdaInelastic* theInelasticModel = new G4LELambdaInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HELambdaInelastic* theHEInelasticModel = new G4HELambdaInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_lambda")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiLambdaInelasticProcess* theInelasticProcess = 
                                new G4AntiLambdaInelasticProcess("inelastic");
      G4LEAntiLambdaInelastic* theInelasticModel =  new G4LEAntiLambdaInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiLambdaInelastic* theHEInelasticModel =  new G4HEAntiLambdaInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4QCaptureAtRest, ordDefault);
    }
    else if (particleName == "sigma+")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4SigmaPlusInelasticProcess* theInelasticProcess = 
                                 new G4SigmaPlusInelasticProcess("inelastic");
      G4LESigmaPlusInelastic* theInelasticModel =  new G4LESigmaPlusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HESigmaPlusInelastic* theHEInelasticModel =  new G4HESigmaPlusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "sigma-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4SigmaMinusInelasticProcess* theInelasticProcess = 
                                 new G4SigmaMinusInelasticProcess("inelastic");
      G4LESigmaMinusInelastic* theInelasticModel =  new G4LESigmaMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HESigmaMinusInelastic* theHEInelasticModel =  new G4HESigmaMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4QCaptureAtRest, ordDefault);
    }
    else if (particleName == "anti_sigma+")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiSigmaPlusInelasticProcess* theInelasticProcess = 
                             new G4AntiSigmaPlusInelasticProcess("inelastic");
      G4LEAntiSigmaPlusInelastic* theInelasticModel =  new G4LEAntiSigmaPlusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiSigmaPlusInelastic* theHEInelasticModel =  new G4HEAntiSigmaPlusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4QCaptureAtRest, ordDefault);
    }
    else if (particleName == "anti_sigma-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiSigmaMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiSigmaMinusInelasticProcess("inelastic");
      G4LEAntiSigmaMinusInelastic* theInelasticModel =  new G4LEAntiSigmaMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiSigmaMinusInelastic* theHEInelasticModel =  new G4HEAntiSigmaMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "xi0")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4XiZeroInelasticProcess* theInelasticProcess = 
                            new G4XiZeroInelasticProcess("inelastic");
      G4LEXiZeroInelastic* theInelasticModel =  new G4LEXiZeroInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEXiZeroInelastic* theHEInelasticModel =  new G4HEXiZeroInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "xi-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4XiMinusInelasticProcess* theInelasticProcess = 
                            new G4XiMinusInelasticProcess("inelastic");
      G4LEXiMinusInelastic* theInelasticModel =  new G4LEXiMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEXiMinusInelastic* theHEInelasticModel =  new G4HEXiMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4QCaptureAtRest, ordDefault);
    }
    else if (particleName == "anti_xi0")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiXiZeroInelasticProcess* theInelasticProcess = 
                            new G4AntiXiZeroInelasticProcess("inelastic");
      G4LEAntiXiZeroInelastic* theInelasticModel =  new G4LEAntiXiZeroInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiXiZeroInelastic* theHEInelasticModel =  new G4HEAntiXiZeroInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "anti_xi-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiXiMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiXiMinusInelasticProcess("inelastic");
      G4LEAntiXiMinusInelastic* theInelasticModel =  new G4LEAntiXiMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiXiMinusInelastic* theHEInelasticModel =  new G4HEAntiXiMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "deuteron")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4DeuteronInelasticProcess* theInelasticProcess = 
                            new G4DeuteronInelasticProcess("inelastic");
      G4LEDeuteronInelastic* theInelasticModel =  new G4LEDeuteronInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "triton")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4TritonInelasticProcess* theInelasticProcess = 
                            new G4TritonInelasticProcess("inelastic");
      G4LETritonInelastic* theInelasticModel =  new G4LETritonInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "alpha")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AlphaInelasticProcess* theInelasticProcess = 
                            new G4AlphaInelasticProcess("inelastic");
      G4LEAlphaInelastic* theInelasticModel =  new G4LEAlphaInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
    else if (particleName == "omega-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4OmegaMinusInelasticProcess* theInelasticProcess = 
                            new G4OmegaMinusInelasticProcess("inelastic");
      G4LEOmegaMinusInelastic* theInelasticModel =  new G4LEOmegaMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEOmegaMinusInelastic* theHEInelasticModel =  new G4HEOmegaMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
      pmanager->AddRestProcess(new G4QCaptureAtRest, ordDefault);
    }
    else if (particleName == "anti_omega-")
    {
      pmanager->AddDiscreteProcess(theElasticProcess);
      G4AntiOmegaMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiOmegaMinusInelasticProcess("inelastic");
      G4LEAntiOmegaMinusInelastic* theInelasticModel =  new G4LEAntiOmegaMinusInelastic;
      theInelasticProcess->RegisterMe(theInelasticModel);
      G4HEAntiOmegaMinusInelastic* theHEInelasticModel =  new G4HEAntiOmegaMinusInelastic;
      theInelasticProcess->RegisterMe(theHEInelasticModel);
      pmanager->AddDiscreteProcess(theInelasticProcess);
    }
  }
}

void PhotInPhysicsList::SetCuts()
{
  if(verboseLevel>0) G4cout<<"PhotInPhysicsList::SetCuts: default cut length : "
                           <<G4BestUnit(defaultCutValue,"Length")<<G4endl;
  // These values are used as the default production thresholds for the world volume.
  SetCutsWithDefault();

  // Production thresholds for detector regions

  G4double fact = 1.; // Multiplicative factor for default cuts
  for(G4int i=0; i< PhotInNumSections; i++)
  { 
    G4Region* reg = G4RegionStore::GetInstance()-> GetRegion(PhotInRegName[i]);
    G4ProductionCuts* cuts = new G4ProductionCuts;
    cuts->SetProductionCut(defaultCutValue*fact);
    reg->SetProductionCuts(cuts);
    fact *= 10.; // @@ Increment the multiplicative factor by order of magnitude
  }
}


