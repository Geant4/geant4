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
// $Id: RE01HadronPhysics.cc,v 1.3 2010-04-07 01:27:53 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "RE01HadronPhysics.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

// processes
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4HadronElasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
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
#include "G4OmegaMinusInelasticProcess.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"
#include "G4PionMinusAbsorptionAtRest.hh"
#include "G4KaonMinusAbsorption.hh"

// cross sections
#include "G4PiNuclearCrossSection.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"

// models
#include "G4LElastic.hh"
#include "G4CascadeInterface.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4LEAntiProtonInelastic.hh"
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
#include "G4LEOmegaMinusInelastic.hh"
#include "G4LEAntiOmegaMinusInelastic.hh"

#include "G4HEAntiProtonInelastic.hh"
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

#include "G4TheoFSGenerator.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4QGSMFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4QGSParticipants.hh"
#include "G4QGSModel.hh"


RE01HadronPhysics::RE01HadronPhysics(const G4String& name)
                          :G4VPhysicsConstructor(name)
{;}

RE01HadronPhysics::~RE01HadronPhysics()
{;}


void RE01HadronPhysics::ConstructParticle()
{
  //  Construct all mesons
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  //  Construct all baryons
  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}


void RE01HadronPhysics::ConstructProcess()
{
  // Hadronic Elastic Process and Model (the same for all hadrons)

  G4HadronElasticProcess* elasticProcess = new G4HadronElasticProcess();
  G4LElastic* elasticModel = new G4LElastic();
  elasticProcess->RegisterMe(elasticModel);

  // Hadronic inelastic models

  // Bertini cascade model: use for p,n,pi+,pi- between 0 and 9.9 GeV
  G4CascadeInterface* bertiniModel = new G4CascadeInterface();
  bertiniModel->SetMaxEnergy(9.9*GeV);

  // Low energy parameterized models : use between 9.5 and 25 GeV
  G4double LEPUpperLimit = 25*GeV;
  G4double LEPpnpiLimit = 9.5*GeV;

  G4LEKaonZeroLInelastic* LEPk0LModel = new G4LEKaonZeroLInelastic();
  LEPk0LModel->SetMaxEnergy(LEPUpperLimit);

  G4LEKaonZeroSInelastic* LEPk0SModel = new G4LEKaonZeroSInelastic();
  LEPk0SModel->SetMaxEnergy(LEPUpperLimit);

  // Quark-Gluon String Model: use for p,n,pi+,pi- between 12 GeV and 100 TeV

  G4TheoFSGenerator* QGSPModel = new G4TheoFSGenerator();
  G4GeneratorPrecompoundInterface* theCascade = 
                                    new G4GeneratorPrecompoundInterface();
  G4ExcitationHandler* exHandler = new G4ExcitationHandler();
  G4PreCompoundModel* preCompound = new G4PreCompoundModel(exHandler);
  theCascade->SetDeExcitation(preCompound);
  QGSPModel->SetTransport(theCascade);
  G4QGSMFragmentation* frag = new G4QGSMFragmentation();
  G4ExcitedStringDecay* stringDecay = new G4ExcitedStringDecay(frag);
  G4QGSModel<G4QGSParticipants>* stringModel = 
                                   new G4QGSModel<G4QGSParticipants>();
  stringModel->SetFragmentationModel(stringDecay);
  QGSPModel->SetHighEnergyGenerator(stringModel);   
  QGSPModel->SetMinEnergy(12*GeV);
  QGSPModel->SetMaxEnergy(100*TeV);

  //  
  G4ProcessManager * pManager = 0;

  ///////////////////
  //               //
  //  pi+ physics  //
  //               //
  ///////////////////

  pManager = G4PionPlus::PionPlus()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);
 
  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4PionPlusInelasticProcess* pipinelProc = new G4PionPlusInelasticProcess();
  G4PiNuclearCrossSection* pion_XC = new G4PiNuclearCrossSection();
  pipinelProc->AddDataSet(pion_XC);
  pipinelProc->RegisterMe(bertiniModel);

  G4LEPionPlusInelastic* LEPpipModel = new G4LEPionPlusInelastic();
  LEPpipModel->SetMinEnergy(LEPpnpiLimit);
  LEPpipModel->SetMaxEnergy(LEPUpperLimit);
  pipinelProc->RegisterMe(LEPpipModel);

  pipinelProc->RegisterMe(QGSPModel);
  pManager->AddDiscreteProcess(pipinelProc);

  ///////////////////
  //               //
  //  pi- physics  //
  //               //
  ///////////////////

  pManager = G4PionMinus::PionMinus()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4PionMinusInelasticProcess* piminelProc = new G4PionMinusInelasticProcess();
  piminelProc->AddDataSet(pion_XC);
  piminelProc->RegisterMe(bertiniModel);

  G4LEPionMinusInelastic* LEPpimModel = new G4LEPionMinusInelastic();
  LEPpimModel->SetMinEnergy(LEPpnpiLimit);
  LEPpimModel->SetMaxEnergy(LEPUpperLimit);
  piminelProc->RegisterMe(LEPpimModel);

  piminelProc->RegisterMe(QGSPModel);
  pManager->AddDiscreteProcess(piminelProc);

  // pi- absorption at rest
  G4PionMinusAbsorptionAtRest* pimAbsorb = new G4PionMinusAbsorptionAtRest();
  pManager->AddRestProcess(pimAbsorb);
   
  ///////////////////
  //               //
  //  K+ physics   //
  //               //
  ///////////////////

  pManager = G4KaonPlus::KaonPlus()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4KaonPlusInelasticProcess* kpinelProc = new G4KaonPlusInelasticProcess();
  G4LEKaonPlusInelastic* LEPkpModel = new G4LEKaonPlusInelastic();
  LEPkpModel->SetMaxEnergy(LEPUpperLimit);
  kpinelProc->RegisterMe(LEPkpModel);
  kpinelProc->RegisterMe(QGSPModel);
  pManager->AddDiscreteProcess(kpinelProc);

  ///////////////////
  //               //
  //  K- physics   //
  //               //
  ///////////////////

  pManager = G4KaonMinus::KaonMinus()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4KaonMinusInelasticProcess* kminelProc = new G4KaonMinusInelasticProcess();
  G4LEKaonMinusInelastic* LEPkmModel = new G4LEKaonMinusInelastic();
  LEPkmModel->SetMaxEnergy(LEPUpperLimit);
  kminelProc->RegisterMe(LEPkmModel);
  kminelProc->RegisterMe(QGSPModel);
  pManager->AddDiscreteProcess(kminelProc);

  // K- absorption at rest
  G4KaonMinusAbsorption* kmAbsorb = new G4KaonMinusAbsorption();
  pManager->AddRestProcess(kmAbsorb);

  ///////////////////
  //               //
  //  K0L physics  //
  //               //
  ///////////////////

  pManager = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4KaonZeroLInelasticProcess* k0LinelProc = new G4KaonZeroLInelasticProcess();
  k0LinelProc->RegisterMe(LEPk0LModel);
  k0LinelProc->RegisterMe(QGSPModel);
  pManager->AddDiscreteProcess(k0LinelProc);

  ///////////////////
  //               //
  //  K0S physics  //
  //               //
  ///////////////////

  pManager = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4KaonZeroSInelasticProcess* k0SinelProc = new G4KaonZeroSInelasticProcess();
  k0SinelProc->RegisterMe(LEPk0SModel);
  k0SinelProc->RegisterMe(QGSPModel);
  pManager->AddDiscreteProcess(k0SinelProc);

  ///////////////////
  //               //
  //    Proton     //
  //               //
  ///////////////////

  pManager = G4Proton::Proton()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4ProtonInelasticProcess* pinelProc = new G4ProtonInelasticProcess();
  G4ProtonInelasticCrossSection* proton_XC = 
                                   new G4ProtonInelasticCrossSection();
  pinelProc->AddDataSet(proton_XC);
  pinelProc->RegisterMe(bertiniModel);

  G4LEProtonInelastic* LEPpModel = new G4LEProtonInelastic();
  LEPpModel->SetMinEnergy(LEPpnpiLimit);
  LEPpModel->SetMaxEnergy(LEPUpperLimit);
  pinelProc->RegisterMe(LEPpModel);

  pinelProc->RegisterMe(QGSPModel);
  pManager->AddDiscreteProcess(pinelProc);

  ///////////////////
  //               //
  //  Anti-Proton  //
  //               //
  ///////////////////

  pManager = G4AntiProton::AntiProton()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiProtonInelasticProcess* apinelProc = 
                                   new G4AntiProtonInelasticProcess();
  G4LEAntiProtonInelastic* LEPapModel = new G4LEAntiProtonInelastic(); 
  apinelProc->RegisterMe(LEPapModel);
  G4HEAntiProtonInelastic* HEPapModel = new G4HEAntiProtonInelastic(); 
  apinelProc->RegisterMe(HEPapModel);
  pManager->AddDiscreteProcess(apinelProc);

  // anti-proton annihilation at rest
  G4AntiProtonAnnihilationAtRest* apAnnihil = 
                                 new G4AntiProtonAnnihilationAtRest();
  pManager->AddRestProcess(apAnnihil);
  
  ///////////////////
  //               //
  //    Neutron    //
  //               //
  ///////////////////

  pManager = G4Neutron::Neutron()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4NeutronInelasticProcess* ninelProc = new G4NeutronInelasticProcess();
  G4NeutronInelasticCrossSection* neutron_XC = 
                                   new G4NeutronInelasticCrossSection();
  ninelProc->AddDataSet(neutron_XC);
  ninelProc->RegisterMe(bertiniModel);

  G4LENeutronInelastic* LEPnModel = new G4LENeutronInelastic();
  LEPnModel->SetMinEnergy(LEPpnpiLimit);
  LEPnModel->SetMaxEnergy(LEPUpperLimit);
  ninelProc->RegisterMe(LEPnModel);

  ninelProc->RegisterMe(QGSPModel);
  pManager->AddDiscreteProcess(ninelProc);

  // neutron-induced fission
  G4HadronFissionProcess* neutronFission = new G4HadronFissionProcess();
  G4LFission* neutronFissionModel = new G4LFission();
  neutronFissionModel->SetMinEnergy(0.);
  neutronFissionModel->SetMaxEnergy(20*TeV);
  neutronFission->RegisterMe(neutronFissionModel);
  pManager->AddDiscreteProcess(neutronFission);

  // neutron capture
  G4HadronCaptureProcess* neutronCapture = new G4HadronCaptureProcess();
  G4LCapture* neutronCaptureModel = new G4LCapture();
  neutronCaptureModel->SetMinEnergy(0.);
  neutronCaptureModel->SetMaxEnergy(20*TeV);
  neutronCapture->RegisterMe(neutronCaptureModel);
  pManager->AddDiscreteProcess(neutronCapture);

  ///////////////////
  //               //
  // Anti-Neutron  //
  //               //
  ///////////////////

  pManager = G4AntiNeutron::AntiNeutron()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiNeutronInelasticProcess* aninelProc = 
                                  new G4AntiNeutronInelasticProcess();
  G4LEAntiNeutronInelastic* LEPanModel = new G4LEAntiNeutronInelastic(); 
  aninelProc->RegisterMe(LEPanModel);
  G4HEAntiNeutronInelastic* HEPanModel = new G4HEAntiNeutronInelastic(); 
  aninelProc->RegisterMe(HEPanModel);
  pManager->AddDiscreteProcess(aninelProc);

  // anti-neutron annihilation at rest
  G4AntiNeutronAnnihilationAtRest* anAnnihil = 
                                 new G4AntiNeutronAnnihilationAtRest();
  pManager->AddRestProcess(anAnnihil);

  ///////////////////
  //               //
  //    Lambda     //
  //               //
  ///////////////////

  pManager = G4Lambda::Lambda()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4LambdaInelasticProcess* linelProc = 
                                  new G4LambdaInelasticProcess();
  G4LELambdaInelastic* LEPlModel = new G4LELambdaInelastic(); 
  linelProc->RegisterMe(LEPlModel);
  G4HELambdaInelastic* HEPlModel = new G4HELambdaInelastic(); 
  linelProc->RegisterMe(HEPlModel);
  pManager->AddDiscreteProcess(linelProc);

  ///////////////////
  //               //
  //  Anti-Lambda  //
  //               //
  ///////////////////

  pManager = G4AntiLambda::AntiLambda()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiLambdaInelasticProcess* alinelProc = 
                                  new G4AntiLambdaInelasticProcess();
  G4LEAntiLambdaInelastic* LEPalModel = new G4LEAntiLambdaInelastic(); 
  alinelProc->RegisterMe(LEPalModel);
  G4HEAntiLambdaInelastic* HEPalModel = new G4HEAntiLambdaInelastic(); 
  alinelProc->RegisterMe(HEPalModel);
  pManager->AddDiscreteProcess(alinelProc);

  ///////////////////
  //               //
  //    Sigma-     //
  //               //
  ///////////////////

  pManager = G4SigmaMinus::SigmaMinus()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4SigmaMinusInelasticProcess* sminelProc = 
                                  new G4SigmaMinusInelasticProcess();
  G4LESigmaMinusInelastic* LEPsmModel = new G4LESigmaMinusInelastic(); 
  sminelProc->RegisterMe(LEPsmModel);
  G4HESigmaMinusInelastic* HEPsmModel = new G4HESigmaMinusInelastic(); 
  sminelProc->RegisterMe(HEPsmModel);
  pManager->AddDiscreteProcess(sminelProc);

  ///////////////////
  //               //
  //  Anti-Sigma-  //
  //               //
  ///////////////////

  pManager = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiSigmaMinusInelasticProcess* asminelProc = 
                                  new G4AntiSigmaMinusInelasticProcess();
  G4LEAntiSigmaMinusInelastic* LEPasmModel = 
                                  new G4LEAntiSigmaMinusInelastic(); 
  asminelProc->RegisterMe(LEPasmModel);
  G4HEAntiSigmaMinusInelastic* HEPasmModel = 
                                  new G4HEAntiSigmaMinusInelastic(); 
  asminelProc->RegisterMe(HEPasmModel);
  pManager->AddDiscreteProcess(asminelProc);

  ///////////////////
  //               //
  //    Sigma+     //
  //               //
  ///////////////////

  pManager = G4SigmaPlus::SigmaPlus()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4SigmaPlusInelasticProcess* spinelProc = new G4SigmaPlusInelasticProcess();
  G4LESigmaPlusInelastic* LEPspModel = new G4LESigmaPlusInelastic(); 
  spinelProc->RegisterMe(LEPspModel);
  G4HESigmaPlusInelastic* HEPspModel = new G4HESigmaPlusInelastic(); 
  spinelProc->RegisterMe(HEPspModel);
  pManager->AddDiscreteProcess(spinelProc);

  ///////////////////
  //               //
  //  Anti-Sigma+  //
  //               //
  ///////////////////

  pManager = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiSigmaPlusInelasticProcess* aspinelProc = 
                                  new G4AntiSigmaPlusInelasticProcess();
  G4LEAntiSigmaPlusInelastic* LEPaspModel = 
                                  new G4LEAntiSigmaPlusInelastic(); 
  aspinelProc->RegisterMe(LEPaspModel);
  G4HEAntiSigmaPlusInelastic* HEPaspModel = 
                                  new G4HEAntiSigmaPlusInelastic(); 
  aspinelProc->RegisterMe(HEPaspModel);
  pManager->AddDiscreteProcess(aspinelProc);

  ///////////////////
  //               //
  //      Xi-      //
  //               //
  ///////////////////

  pManager = G4XiMinus::XiMinus()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4XiMinusInelasticProcess* xminelProc = new G4XiMinusInelasticProcess();
  G4LEXiMinusInelastic* LEPxmModel = new G4LEXiMinusInelastic(); 
  xminelProc->RegisterMe(LEPxmModel);
  G4HEXiMinusInelastic* HEPxmModel = new G4HEXiMinusInelastic(); 
  xminelProc->RegisterMe(HEPxmModel);
  pManager->AddDiscreteProcess(xminelProc);

  ///////////////////
  //               //
  //   Anti-Xi-    //
  //               //
  ///////////////////

  pManager = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiXiMinusInelasticProcess* axminelProc = 
                                  new G4AntiXiMinusInelasticProcess();
  G4LEAntiXiMinusInelastic* LEPaxmModel = new G4LEAntiXiMinusInelastic(); 
  axminelProc->RegisterMe(LEPaxmModel);
  G4HEAntiXiMinusInelastic* HEPaxmModel = new G4HEAntiXiMinusInelastic(); 
  axminelProc->RegisterMe(HEPaxmModel);
  pManager->AddDiscreteProcess(axminelProc);

  ///////////////////
  //               //
  //      Xi0      //
  //               //
  ///////////////////

  pManager = G4XiZero::XiZero()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4XiZeroInelasticProcess* x0inelProc = new G4XiZeroInelasticProcess();
  G4LEXiZeroInelastic* LEPx0Model = new G4LEXiZeroInelastic(); 
  x0inelProc->RegisterMe(LEPx0Model);
  G4HEXiZeroInelastic* HEPx0Model = new G4HEXiZeroInelastic(); 
  x0inelProc->RegisterMe(HEPx0Model);
  pManager->AddDiscreteProcess(x0inelProc);

  ///////////////////
  //               //
  //   Anti-Xi0    //
  //               //
  ///////////////////

  pManager = G4AntiXiZero::AntiXiZero()->GetProcessManager();

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiXiZeroInelasticProcess* ax0inelProc = 
                                new G4AntiXiZeroInelasticProcess();
  G4LEAntiXiZeroInelastic* LEPax0Model = new G4LEAntiXiZeroInelastic(); 
  ax0inelProc->RegisterMe(LEPax0Model);
  G4HEAntiXiZeroInelastic* HEPax0Model = new G4HEAntiXiZeroInelastic(); 
  ax0inelProc->RegisterMe(HEPax0Model);
  pManager->AddDiscreteProcess(ax0inelProc);

  ///////////////////
  //               //
  //    Omega-     //
  //               //
  ///////////////////

  pManager = G4OmegaMinus::OmegaMinus()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4OmegaMinusInelasticProcess* ominelProc = 
                                      new G4OmegaMinusInelasticProcess();
  G4LEOmegaMinusInelastic* LEPomModel = new G4LEOmegaMinusInelastic(); 
  ominelProc->RegisterMe(LEPomModel);
  G4HEOmegaMinusInelastic* HEPomModel = new G4HEOmegaMinusInelastic(); 
  ominelProc->RegisterMe(HEPomModel);
  pManager->AddDiscreteProcess(ominelProc);

  ///////////////////
  //               //
  //  Anti-Omega-  //
  //               //
  ///////////////////

  pManager = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();

  // EM processes
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

  // hadron elastic
  pManager->AddDiscreteProcess(elasticProcess);

  // hadron inelastic
  G4AntiOmegaMinusInelasticProcess* aominelProc = 
                                      new G4AntiOmegaMinusInelasticProcess();
  G4LEAntiOmegaMinusInelastic* LEPaomModel = 
                                      new G4LEAntiOmegaMinusInelastic(); 
  aominelProc->RegisterMe(LEPaomModel);
  G4HEAntiOmegaMinusInelastic* HEPaomModel = 
                                      new G4HEAntiOmegaMinusInelastic(); 
  aominelProc->RegisterMe(HEPaomModel);
  pManager->AddDiscreteProcess(aominelProc);

}
