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
/// \file runAndEvent/RE02/src/RE02HadronPhysics.cc
/// \brief Implementation of the RE02HadronPhysics class
//
// $Id$
// --------------------------------------------------------------
//
// 09-Oct-2003 Hadron Physics List with Parameterization Model by T. Koi
// 12-Oct-2003 Bug Fixed (KaonMinus) by T. Koi
// 16-Nov-2005 Binary Cascade for Protons. by T. Aso
//             Proton :  BinaryCascade < 6 GeV, 4 GeV < LE Model, 
//             Neutron:  BinaryCascade < 6 GeV, 4 GeV < LE Model, 
//

#include "RE02HadronPhysics.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"    
#include "G4ios.hh"
#include <iomanip>   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE02HadronPhysics::RE02HadronPhysics(const G4String& name)
                   :  G4VPhysicsConstructor(name)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE02HadronPhysics::~RE02HadronPhysics()
{
}

#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4ProcessManager.hh"


void RE02HadronPhysics::ConstructProcess()
{
  G4ProcessManager* pManager = 0;

  // Set up models
  G4HadronElastic* elasticModel = new G4HadronElastic;
  G4CascadeInterface* bertini = new G4CascadeInterface;
  bertini->SetMinEnergy(0);
  bertini->SetMaxEnergy(9*GeV);

  // FTF for high energy
  G4TheoFSGenerator* theTheoModel = new G4TheoFSGenerator;
  G4FTFModel* theFTFModel = new G4FTFModel;
  G4ExcitedStringDecay* stringDecay =
                       new G4ExcitedStringDecay(new G4LundStringFragmentation);
  theFTFModel->SetFragmentationModel(stringDecay);

  G4GeneratorPrecompoundInterface* precoCascade =
                                           new G4GeneratorPrecompoundInterface;
  G4PreCompoundModel* preEquilib =
                               new G4PreCompoundModel(new G4ExcitationHandler);
  precoCascade->SetDeExcitation(preEquilib);
  theTheoModel->SetTransport(precoCascade);
  theTheoModel->SetHighEnergyGenerator(theFTFModel);
  theTheoModel->SetMinEnergy(6*GeV);
  theTheoModel->SetMaxEnergy(20*TeV);

  // FTF for anti-baryons at all energies
  G4TheoFSGenerator* antiBHighEnergyModel = new G4TheoFSGenerator("ANTI-FTFP");
  G4FTFModel* antiBStringModel = new G4FTFModel;
  G4ExcitedStringDecay* stringDecayforAntiB =
                       new G4ExcitedStringDecay(new G4LundStringFragmentation);
  antiBStringModel->SetFragmentationModel(stringDecayforAntiB);

  G4GeneratorPrecompoundInterface* antiBCascade =
                                           new G4GeneratorPrecompoundInterface;
  G4PreCompoundModel* preEquilibforAntiB =
                               new G4PreCompoundModel(new G4ExcitationHandler);
  antiBCascade->SetDeExcitation(preEquilibforAntiB);
  antiBHighEnergyModel->SetTransport(antiBCascade);
  antiBHighEnergyModel->SetHighEnergyGenerator(antiBStringModel);
  antiBHighEnergyModel->SetMinEnergy(0.0);
  antiBHighEnergyModel->SetMaxEnergy(20*TeV);

   
   // Pi+ Physics
   pManager = G4PionPlus::PionPlus()->GetProcessManager();

   // add processes
   G4HadronElasticProcess* theppElasticProcess = new G4HadronElasticProcess();
   theppElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(theppElasticProcess);
 
   G4PionPlusInelasticProcess* thePionPlusInelasticProcess =
                                             new G4PionPlusInelasticProcess(); 
   thePionPlusInelasticProcess->RegisterMe(bertini);
   thePionPlusInelasticProcess->RegisterMe(theTheoModel);
   pManager->AddDiscreteProcess(thePionPlusInelasticProcess);

   G4VProcess* theppMultipleScattering = new G4hMultipleScattering();
   G4VProcess* theppIonisation         = new G4hIonisation();
   // 
   pManager->AddProcess(theppIonisation);
   pManager->AddProcess(theppMultipleScattering);
   // 
   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(theppMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(theppIonisation,         idxAlongStep,2);
   // 
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(theppMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(theppIonisation,         idxPostStep,2);

   // Pi- Physics
   pManager = G4PionMinus::PionMinus()->GetProcessManager();

   G4HadronElasticProcess* thepmElasticProcess = new G4HadronElasticProcess();
   thepmElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(thepmElasticProcess);

   G4PionMinusInelasticProcess* thePionMinusInelasticProcess =
                                            new G4PionMinusInelasticProcess(); 
   thePionMinusInelasticProcess->RegisterMe(bertini);
   thePionMinusInelasticProcess->RegisterMe(theTheoModel);
   pManager->AddDiscreteProcess(thePionMinusInelasticProcess);

   G4VProcess* thepmMultipleScattering = new G4hMultipleScattering();
   G4VProcess* thepmIonisation         = new G4hIonisation();
   // 
   // add processes
   pManager->AddProcess(thepmIonisation);
   pManager->AddProcess(thepmMultipleScattering);
   // 
   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(thepmMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(thepmIonisation,         idxAlongStep,2);
   // 
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(thepmMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(thepmIonisation,         idxPostStep,2);

   // K+ Physics
   pManager = G4KaonPlus::KaonPlus()->GetProcessManager();

   G4HadronElasticProcess* thekpElasticProcess = new G4HadronElasticProcess();
   thekpElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(thekpElasticProcess);

   G4KaonPlusInelasticProcess* theKaonPlusInelasticProcess =
                                             new G4KaonPlusInelasticProcess(); 
   theKaonPlusInelasticProcess->RegisterMe(bertini);
   theKaonPlusInelasticProcess->RegisterMe(theTheoModel);
   pManager->AddDiscreteProcess(theKaonPlusInelasticProcess);

   G4VProcess* thekpMultipleScattering = new G4hMultipleScattering();
   G4VProcess* thekpIonisation         = new G4hIonisation();
   // 
   // add processes
   pManager->AddProcess(thekpIonisation);
   pManager->AddProcess(thekpMultipleScattering);
   // 
   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(thekpMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(thekpIonisation,         idxAlongStep,2);
   // 
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(thekpMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(thekpIonisation,         idxPostStep,2);

   // K- Physics
   pManager = G4KaonMinus::KaonMinus()->GetProcessManager();

   // add processes
   G4HadronElasticProcess* thekmElasticProcess = new G4HadronElasticProcess();
   thekmElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(thekmElasticProcess);

   G4KaonMinusInelasticProcess* theKaonMinusInelasticProcess =
                                            new G4KaonMinusInelasticProcess(); 
   theKaonMinusInelasticProcess->RegisterMe(bertini);
   theKaonMinusInelasticProcess->RegisterMe(theTheoModel);
   pManager->AddDiscreteProcess(theKaonMinusInelasticProcess);

   G4VProcess* thekmMultipleScattering = new G4hMultipleScattering();
   G4VProcess* thekmIonisation         = new G4hIonisation();

   pManager->AddProcess(thekmIonisation);
   pManager->AddProcess(thekmMultipleScattering);

   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(thekmMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(thekmIonisation,         idxAlongStep,2);
   // 
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(thekmMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(thekmIonisation,         idxPostStep,2);

   // Kaon0L Physics
   pManager = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();

   G4HadronElasticProcess* thek0lElasticProcess = new G4HadronElasticProcess();
   thek0lElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(thek0lElasticProcess);

   G4KaonZeroLInelasticProcess* theKaonZeroLInelasticProcess =
                                             new G4KaonZeroLInelasticProcess(); 
   theKaonZeroLInelasticProcess->RegisterMe(bertini);
   theKaonZeroLInelasticProcess->RegisterMe(theTheoModel);
   pManager->AddDiscreteProcess(theKaonZeroLInelasticProcess);

   // Kaon0S Physics
   pManager = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();

   G4HadronElasticProcess* thek0sElasticProcess = new G4HadronElasticProcess();
   thek0sElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(thek0sElasticProcess);

   G4KaonZeroSInelasticProcess* theKaonZeroSInelasticProcess =
                                             new G4KaonZeroSInelasticProcess(); 
   theKaonZeroSInelasticProcess->RegisterMe(bertini);
   theKaonZeroSInelasticProcess->RegisterMe(theTheoModel);
   pManager->AddDiscreteProcess(theKaonZeroSInelasticProcess);

   // Proton Physics
   pManager = G4Proton::Proton()->GetProcessManager();

   // add process
   G4HadronElasticProcess* thepElasticProcess = new G4HadronElasticProcess();
   thepElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(thepElasticProcess);

   G4ProtonInelasticProcess* theProtonInelasticProcess =
     new G4ProtonInelasticProcess(); 

   G4BinaryCascade* theProtonBCModel = new G4BinaryCascade();
   theProtonBCModel->SetMaxEnergy(7.*GeV);
   theProtonInelasticProcess->RegisterMe(theProtonBCModel);
   theProtonInelasticProcess->RegisterMe(theTheoModel);
   pManager->AddDiscreteProcess(theProtonInelasticProcess);

   G4VProcess* thepMultipleScattering = new G4hMultipleScattering();
   G4VProcess* thepIonisation         = new G4hIonisation();

   pManager->AddProcess(thepIonisation);
   pManager->AddProcess(thepMultipleScattering);

   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(thepMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(thepIonisation,         idxAlongStep,2);
   // 
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(thepMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(thepIonisation,         idxPostStep,2);

   // anti-proton Physics
   pManager = G4AntiProton::AntiProton()->GetProcessManager();

   // add process
   G4HadronElasticProcess* theapElasticProcess = new G4HadronElasticProcess();
   theapElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(theapElasticProcess);

   G4AntiProtonInelasticProcess* theAntiProtonInelasticProcess =
                                           new G4AntiProtonInelasticProcess(); 
   theAntiProtonInelasticProcess->RegisterMe(antiBHighEnergyModel);
   pManager->AddDiscreteProcess(theAntiProtonInelasticProcess);

   G4AntiProtonAbsorptionFritiof* theAntiProtonAnnihilation =
      new G4AntiProtonAbsorptionFritiof();
   pManager->AddRestProcess(theAntiProtonAnnihilation);

   G4VProcess* theapMultipleScattering = new G4hMultipleScattering();
   G4VProcess* theapIonisation         = new G4hIonisation();

   pManager->AddProcess(theapIonisation);
   pManager->AddProcess(theapMultipleScattering);

   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(theapMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(theapIonisation,         idxAlongStep,2);
   // 
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(theapMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(theapIonisation,         idxPostStep,2);

   // neutron Physics
   pManager = G4Neutron::Neutron()->GetProcessManager();

   // add process
   G4HadronElasticProcess* thenElasticProcess = new G4HadronElasticProcess();
   thenElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(thenElasticProcess);

   G4NeutronInelasticProcess* theNeutronInelasticProcess =
     new G4NeutronInelasticProcess(); 

   G4BinaryCascade* theNeutronBCModel = new G4BinaryCascade();
   theNeutronBCModel->SetMaxEnergy(7.*GeV);
   theNeutronInelasticProcess->RegisterMe(theNeutronBCModel);
   theNeutronInelasticProcess->RegisterMe(theTheoModel);
   pManager->AddDiscreteProcess(theNeutronInelasticProcess);

   G4HadronFissionProcess* thenFission = new G4HadronFissionProcess();
   G4LFission* thenFissionModel = new G4LFission();
   thenFission->RegisterMe(thenFissionModel);
   pManager->AddDiscreteProcess(thenFission);

   G4HadronCaptureProcess* thenCapture = new G4HadronCaptureProcess();
   G4NeutronRadCapture* thenCaptureModel = new G4NeutronRadCapture();
   thenCapture->RegisterMe(thenCaptureModel);
   G4NeutronCaptureXS* thencapXS = new G4NeutronCaptureXS;
   thenCapture->AddDataSet(thencapXS);
   pManager->AddDiscreteProcess(thenCapture);

   // anti-neutron Physics
   pManager = G4AntiNeutron::AntiNeutron()->GetProcessManager();

   // add process
   G4HadronElasticProcess* theanElasticProcess = new G4HadronElasticProcess();
   theanElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(theanElasticProcess);

   G4AntiNeutronInelasticProcess* theAntiNeutronInelasticProcess =
                                          new G4AntiNeutronInelasticProcess(); 
   theAntiNeutronInelasticProcess->RegisterMe(antiBHighEnergyModel);
   pManager->AddDiscreteProcess(theAntiNeutronInelasticProcess);

   G4AntiNeutronAnnihilationAtRest* theAntiNeutronAnnihilation =
     new G4AntiNeutronAnnihilationAtRest();
   pManager->AddRestProcess(theAntiNeutronAnnihilation);

   // Lambda Physics
   pManager = G4Lambda::Lambda()->GetProcessManager();

   // add process
   G4HadronElasticProcess* thel0ElasticProcess = new G4HadronElasticProcess();
   thel0ElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(thel0ElasticProcess);

   G4LambdaInelasticProcess* theLambdaInelasticProcess =
                                               new G4LambdaInelasticProcess(); 
   theLambdaInelasticProcess->RegisterMe(bertini);
   theLambdaInelasticProcess->RegisterMe(theTheoModel);
   pManager->AddDiscreteProcess(theLambdaInelasticProcess);

   // Anti-Lambda Physics
   pManager = G4AntiLambda::AntiLambda()->GetProcessManager();

   // add process
   G4HadronElasticProcess* theal0ElasticProcess = new G4HadronElasticProcess();
   theal0ElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(theal0ElasticProcess);

   G4AntiLambdaInelasticProcess* theAntiLambdaInelasticProcess =
                                            new G4AntiLambdaInelasticProcess(); 
   theAntiLambdaInelasticProcess->RegisterMe(antiBHighEnergyModel);
   pManager->AddDiscreteProcess(theAntiLambdaInelasticProcess);

   // Sigma+ Physics
   pManager = G4SigmaPlus::SigmaPlus()->GetProcessManager();

   // add process
   G4HadronElasticProcess* thespElasticProcess = new G4HadronElasticProcess();
   thespElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(thespElasticProcess);

   G4SigmaPlusInelasticProcess* theSigmaPlusInelasticProcess =
                                            new G4SigmaPlusInelasticProcess(); 
   theSigmaPlusInelasticProcess->RegisterMe(bertini);
   theSigmaPlusInelasticProcess->RegisterMe(theTheoModel);
   pManager->AddDiscreteProcess(theSigmaPlusInelasticProcess);

   G4VProcess* thespMultipleScattering = new G4hMultipleScattering();
   G4VProcess* thespIonisation         = new G4hIonisation();

   pManager->AddProcess(thespIonisation);
   pManager->AddProcess(thespMultipleScattering);

   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(thespMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(thespIonisation,         idxAlongStep,2);
   // 
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(thespMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(thespIonisation,         idxPostStep,2);

   // anti-Sigma+ Physics
   pManager = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();

   // add process
   G4HadronElasticProcess* theaspElasticProcess = new G4HadronElasticProcess();
   theaspElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(theaspElasticProcess);

   G4AntiSigmaPlusInelasticProcess* theAntiSigmaPlusInelasticProcess =
                                         new G4AntiSigmaPlusInelasticProcess(); 
   theAntiSigmaPlusInelasticProcess->RegisterMe(antiBHighEnergyModel);
   pManager->AddDiscreteProcess(theAntiSigmaPlusInelasticProcess);

   G4VProcess* theaspMultipleScattering = new G4hMultipleScattering();
   G4VProcess* theaspIonisation         = new G4hIonisation();

   pManager->AddProcess(theaspIonisation);
   pManager->AddProcess(theaspMultipleScattering);

   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(theaspMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(theaspIonisation,         idxAlongStep,2);
   // 
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(theaspMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(theaspIonisation,         idxPostStep,2);

   // Sigma- Physics
   pManager = G4SigmaMinus::SigmaMinus()->GetProcessManager();

   // add process
   G4HadronElasticProcess* thesmElasticProcess = new G4HadronElasticProcess();
   thesmElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(thesmElasticProcess);

   G4SigmaMinusInelasticProcess* theSigmaMinusInelasticProcess =
                                           new G4SigmaMinusInelasticProcess(); 
   theSigmaMinusInelasticProcess->RegisterMe(bertini);
   theSigmaMinusInelasticProcess->RegisterMe(theTheoModel);
   pManager->AddDiscreteProcess(theSigmaMinusInelasticProcess);

   G4VProcess* thesmMultipleScattering = new G4hMultipleScattering();
   G4VProcess* thesmIonisation         = new G4hIonisation();

   pManager->AddProcess(thesmIonisation);
   pManager->AddProcess(thesmMultipleScattering);

   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(thesmMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(thesmIonisation,         idxAlongStep,2);
   // 
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(thesmMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(thesmIonisation,         idxPostStep,2);

   // anti-Sigma- Physics
   pManager = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();

   // add process
   G4HadronElasticProcess* theasmElasticProcess = new G4HadronElasticProcess();
   theasmElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(theasmElasticProcess);

   G4AntiSigmaMinusInelasticProcess* theAntiSigmaMinusInelasticProcess =
     new G4AntiSigmaMinusInelasticProcess(); 
   theAntiSigmaMinusInelasticProcess->RegisterMe(antiBHighEnergyModel);
   pManager->AddDiscreteProcess(theAntiSigmaMinusInelasticProcess);

   G4VProcess* theasmMultipleScattering = new G4hMultipleScattering();
   G4VProcess* theasmIonisation         = new G4hIonisation();

   pManager->AddProcess(theasmIonisation);
   pManager->AddProcess(theasmMultipleScattering);

   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(theasmMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(theasmIonisation,         idxAlongStep,2);
   // 
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(theasmMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(theasmIonisation,         idxPostStep,2);

   // Xi0 Physics
   pManager = G4XiZero::XiZero()->GetProcessManager();

   // add process
   G4HadronElasticProcess* thex0ElasticProcess = new G4HadronElasticProcess();
   thex0ElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(thex0ElasticProcess);

   G4XiZeroInelasticProcess* theXiZeroInelasticProcess =
                                               new G4XiZeroInelasticProcess(); 
   theXiZeroInelasticProcess->RegisterMe(bertini);
   theXiZeroInelasticProcess->RegisterMe(theTheoModel);
   pManager->AddDiscreteProcess(theXiZeroInelasticProcess);

   // Anti-Xi0 Physics
   pManager = G4AntiXiZero::AntiXiZero()->GetProcessManager();

   // add process
   G4HadronElasticProcess* theax0ElasticProcess = new G4HadronElasticProcess();
   theax0ElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(theax0ElasticProcess);

   G4AntiXiZeroInelasticProcess* theAntiXiZeroInelasticProcess =
     new G4AntiXiZeroInelasticProcess(); 
   theAntiXiZeroInelasticProcess->RegisterMe(antiBHighEnergyModel);
   pManager->AddDiscreteProcess(theAntiXiZeroInelasticProcess);

   // Xi- Physics
   pManager = G4XiMinus::XiMinus()->GetProcessManager();

   // add process
   G4HadronElasticProcess* thexmElasticProcess = new G4HadronElasticProcess();
   thexmElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(thexmElasticProcess);

   G4XiMinusInelasticProcess* theXiMinusInelasticProcess =
                                              new G4XiMinusInelasticProcess(); 
   theXiMinusInelasticProcess->RegisterMe(bertini);
   theXiMinusInelasticProcess->RegisterMe(theTheoModel);
   pManager->AddDiscreteProcess(theXiMinusInelasticProcess);

   G4VProcess* thexmMultipleScattering = new G4hMultipleScattering();
   G4VProcess* thexmIonisation         = new G4hIonisation();

   pManager->AddProcess(thexmIonisation);
   pManager->AddProcess(thexmMultipleScattering);

   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(thexmMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(thexmIonisation,         idxAlongStep,2);
   // 
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(thexmMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(thexmIonisation,         idxPostStep,2);

   // anti-Xi- Physics
   pManager = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();

   // add process
   G4HadronElasticProcess* theaxmElasticProcess = new G4HadronElasticProcess();
   theaxmElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(theaxmElasticProcess);

   G4AntiXiMinusInelasticProcess* theAntiXiMinusInelasticProcess =
                                           new G4AntiXiMinusInelasticProcess(); 
   theAntiXiMinusInelasticProcess->RegisterMe(antiBHighEnergyModel);
   pManager->AddDiscreteProcess(theAntiXiMinusInelasticProcess);

   G4VProcess* theaxmMultipleScattering = new G4hMultipleScattering();
   G4VProcess* theaxmIonisation         = new G4hIonisation();

   pManager->AddProcess(theaxmIonisation);
   pManager->AddProcess(theaxmMultipleScattering);

   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(theaxmMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(theaxmIonisation,         idxAlongStep,2);
   // 
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(theaxmMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(theaxmIonisation,         idxPostStep,2);

   // Omega- Physics
   pManager = G4OmegaMinus::OmegaMinus()->GetProcessManager();

   // add process
   G4HadronElasticProcess* theomElasticProcess = new G4HadronElasticProcess();
   theomElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(theomElasticProcess);

   G4OmegaMinusInelasticProcess* theOmegaMinusInelasticProcess =
                                           new G4OmegaMinusInelasticProcess(); 
   theOmegaMinusInelasticProcess->RegisterMe(bertini);
   theOmegaMinusInelasticProcess->RegisterMe(theTheoModel);
   pManager->AddDiscreteProcess(theOmegaMinusInelasticProcess);

   G4VProcess* theomMultipleScattering = new G4hMultipleScattering();
   G4VProcess* theomIonisation         = new G4hIonisation();

   pManager->AddProcess(theomIonisation);
   pManager->AddProcess(theomMultipleScattering);

   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(theomMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(theomIonisation,         idxAlongStep,2);
   // 
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(theomMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(theomIonisation,         idxPostStep,2);

   // anti-Omega- Physics
   pManager = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();

   // add process
   G4HadronElasticProcess* theaomElasticProcess = new G4HadronElasticProcess();
   theaomElasticProcess->RegisterMe(elasticModel);
   pManager->AddDiscreteProcess(theaomElasticProcess);

   G4AntiOmegaMinusInelasticProcess* theAntiOmegaMinusInelasticProcess =
                                        new G4AntiOmegaMinusInelasticProcess(); 
   theAntiOmegaMinusInelasticProcess->RegisterMe(antiBHighEnergyModel);
   pManager->AddDiscreteProcess(theAntiOmegaMinusInelasticProcess);

   G4VProcess* theaomMultipleScattering = new G4hMultipleScattering();
   G4VProcess* theaomIonisation         = new G4hIonisation();

   pManager->AddProcess(theaomIonisation);
   pManager->AddProcess(theaomMultipleScattering);

   // set ordering for AlongStepDoIt
   pManager->SetProcessOrdering(theaomMultipleScattering, idxAlongStep,1);
   pManager->SetProcessOrdering(theaomIonisation,         idxAlongStep,2);
   // 
   // set ordering for PostStepDoIt
   pManager->SetProcessOrdering(theaomMultipleScattering, idxPostStep,1);
   pManager->SetProcessOrdering(theaomIonisation,         idxPostStep,2);

}
