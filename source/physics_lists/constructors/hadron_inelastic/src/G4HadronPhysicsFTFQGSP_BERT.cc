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
// $Id: $
//
//---------------------------------------------------------------------------
// Author: Alberto Ribon
// Date:   October 2017
//
// Hadron physics for the new, experimental physics list FTFQGSP_BERT,
// with QGS fragmentation of strings, instead of the Lund string
// fragmentation. Note that the string excitation is still done with FTF,
// exactly as for FTFP_BERT.
// Given that it is an experimental, and perhaps temporary, new type of
// hadron physics, corresponding builders are not created and everything
// is implemented directly in this class.
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsFTFQGSP_BERT.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4ProcessManager.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4CrossSectionInelastic.hh"

#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTFQGSP_BERT);

G4HadronPhysicsFTFQGSP_BERT::G4HadronPhysicsFTFQGSP_BERT(G4int)
    : G4VPhysicsConstructor("hInelastic FTFQGSP_BERT")
    , theNeutronCaptureModel(0)
    , thePreEquilib(0)
    , theCascade(0)
    , theStringModel(0)
    , theStringDecay(0)
    , theQGSMFragmentation(0)
    , theHandler(0)
    , theModel1(0)
    , theModel2(0)
    , theModel3(0)
    , theBertini1(0)
    , theBertini2(0)
    , theNeutronCaptureProcess(0)
    , theNeutronInelastic(0)
    , theProtonInelastic(0)
    , thePionMinusInelastic(0)
    , thePionPlusInelastic(0)
    , theKaonMinusInelastic(0)
    , theKaonPlusInelastic(0)
    , theKaonZeroLInelastic(0)
    , theKaonZeroSInelastic(0)
    , theLambdaInelastic(0)
    , theAntiLambdaInelastic(0)
    , theSigmaMinusInelastic(0)
    , theAntiSigmaMinusInelastic(0)
    , theSigmaPlusInelastic(0)
    , theAntiSigmaPlusInelastic(0)
    , theXiZeroInelastic(0)
    , theAntiXiZeroInelastic(0)
    , theXiMinusInelastic(0)
    , theAntiXiMinusInelastic(0)
    , theOmegaMinusInelastic(0)
    , theAntiOmegaMinusInelastic(0)
    , theAntiProtonInelastic(0)
    , theAntiNeutronInelastic(0)
    , theAntiDeuteronInelastic(0)
    , theAntiTritonInelastic(0)
    , theAntiHe3Inelastic(0)
    , theAntiAlphaInelastic(0)
    , thePiXS(0)
    , theKaonXS(0)
    , theChipsHyperonInelasticXS(0)
    , theAntiNucleonXS(0)
    , theNeutronInelasticXS(0)
    , theNeutronCaptureXS(0)
{}


G4HadronPhysicsFTFQGSP_BERT::G4HadronPhysicsFTFQGSP_BERT(const G4String& name, G4bool /* quasiElastic */)
    : G4VPhysicsConstructor(name) 
    , theNeutronCaptureModel(0)
    , thePreEquilib(0)
    , theCascade(0)
    , theStringModel(0)
    , theStringDecay(0)
    , theQGSMFragmentation(0)
    , theHandler(0)
    , theModel1(0)
    , theModel2(0)
    , theModel3(0)
    , theBertini1(0)
    , theBertini2(0)
    , theNeutronCaptureProcess(0)
    , theNeutronInelastic(0)
    , theProtonInelastic(0)
    , thePionMinusInelastic(0)
    , thePionPlusInelastic(0)
    , theKaonMinusInelastic(0)
    , theKaonPlusInelastic(0)
    , theKaonZeroLInelastic(0)
    , theKaonZeroSInelastic(0)
    , theLambdaInelastic(0)
    , theAntiLambdaInelastic(0)
    , theSigmaMinusInelastic(0)
    , theAntiSigmaMinusInelastic(0)
    , theSigmaPlusInelastic(0)
    , theAntiSigmaPlusInelastic(0)
    , theXiZeroInelastic(0)
    , theAntiXiZeroInelastic(0)
    , theXiMinusInelastic(0)
    , theAntiXiMinusInelastic(0)
    , theOmegaMinusInelastic(0)
    , theAntiOmegaMinusInelastic(0)
    , theAntiProtonInelastic(0)
    , theAntiNeutronInelastic(0)
    , theAntiDeuteronInelastic(0)
    , theAntiTritonInelastic(0)
    , theAntiHe3Inelastic(0)
    , theAntiAlphaInelastic(0)
    , thePiXS(0)
    , theKaonXS(0)
    , theChipsHyperonInelasticXS(0)
    , theAntiNucleonXS(0)
    , theNeutronInelasticXS(0)
    , theNeutronCaptureXS(0)
{}


void G4HadronPhysicsFTFQGSP_BERT::CreateModels()
{

  G4double minFTFP =  3.0 * GeV;
  G4double maxBERT = 12.0 * GeV;
  G4cout << " FTFQGSP_BERT : similar to FTFP_BERT but with" << G4endl
         << " QGS string fragmentation (instead of Lund string fragmentation)." << G4endl;

  theStringModel = new G4FTFModel;

  theStringDecay = new G4ExcitedStringDecay( theQGSMFragmentation = new G4QGSMFragmentation );
  theStringModel->SetFragmentationModel( theStringDecay );
  thePreEquilib = new G4PreCompoundModel( theHandler = new G4ExcitationHandler );
  theCascade = new G4GeneratorPrecompoundInterface( thePreEquilib );

  // FTF for neutrons, protons, pions, and kaons
  theModel1 = new G4TheoFSGenerator( "FTFP" );
  theModel1->SetMinEnergy( minFTFP );
  theModel1->SetMaxEnergy( 100.0*TeV );
  theModel1->SetTransport( theCascade );
  theModel1->SetHighEnergyGenerator( theStringModel );
 
  // BERT for neutrons, protons, pions, and kaons
  theBertini1 = new G4CascadeInterface;
  theBertini1->SetMinEnergy( 0.0*GeV );
  theBertini1->SetMaxEnergy( maxBERT );

  // FTF for hyperons
  theModel2 = new G4TheoFSGenerator( "FTFP" );
  theModel2->SetMinEnergy( 2.0*GeV );
  theModel2->SetMaxEnergy( 100.0*TeV );
  theModel2->SetTransport( theCascade );
  theModel2->SetHighEnergyGenerator( theStringModel );
  
  // BERT for hyperons
  theBertini2 = new G4CascadeInterface;
  theBertini2->SetMinEnergy( 0.0*GeV );
  theBertini2->SetMaxEnergy( 6.0*GeV );

  // FTF for Antibaryons  
  theModel3 = new G4TheoFSGenerator( "FTFP" );
  theModel3->SetMinEnergy( 0.0*GeV );
  theModel3->SetMaxEnergy( 100.0*TeV );
  theModel3->SetTransport( theCascade );
  theModel3->SetHighEnergyGenerator( theStringModel );

  // Neutron Capture
  theNeutronCaptureModel = new G4NeutronRadCapture;
  theNeutronCaptureModel->SetMinEnergy( 0.0 );
  theNeutronCaptureModel->SetMaxEnergy( 100.0*TeV );

  // Cross sections
  thePiXS = new G4CrossSectionPairGG( new G4PiNuclearCrossSection, 91*GeV );
  theAntiNucleonXS = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
  theKaonXS = new G4CrossSectionInelastic( new G4ComponentGGHadronNucleusXsc );
  theChipsHyperonInelasticXS = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet( G4ChipsHyperonInelasticXS::Default_Name() );
  theNeutronInelasticXS      = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet( G4NeutronInelasticXS::Default_Name() );
  theNeutronCaptureXS        = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet( G4NeutronCaptureXS::Default_Name() );
}


G4HadronPhysicsFTFQGSP_BERT::~G4HadronPhysicsFTFQGSP_BERT()
{
  delete theStringDecay;
  delete theStringModel;
  delete thePreEquilib;
  delete theCascade;
  delete theQGSMFragmentation;
}


void G4HadronPhysicsFTFQGSP_BERT::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}


#include "G4ProcessManager.hh"
void G4HadronPhysicsFTFQGSP_BERT::ConstructProcess()
{
  CreateModels();

  G4ProcessManager * aProcMan = 0;

  theNeutronInelastic = new G4NeutronInelasticProcess();
  theNeutronInelastic->RegisterMe( theModel1 );
  theNeutronInelastic->RegisterMe( theBertini1 );
  theNeutronInelastic->AddDataSet( new G4BGGNucleonInelasticXS( G4Neutron::Neutron() ) );
  theNeutronInelastic->AddDataSet( theNeutronInelasticXS );
  aProcMan = G4Neutron::Neutron()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theNeutronInelastic );
  theNeutronCaptureProcess = new G4HadronCaptureProcess();
  theNeutronCaptureProcess->RegisterMe( theNeutronCaptureModel );
  theNeutronCaptureProcess->AddDataSet( theNeutronCaptureXS );
  aProcMan->AddDiscreteProcess( theNeutronCaptureProcess );

  theProtonInelastic = new G4ProtonInelasticProcess();
  theProtonInelastic->RegisterMe( theModel1 );
  theProtonInelastic->RegisterMe( theBertini1 );
  theProtonInelastic->AddDataSet( new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
  aProcMan = G4Proton::Proton()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theProtonInelastic );

  thePionMinusInelastic = new G4PionMinusInelasticProcess();
  thePionMinusInelastic->RegisterMe( theModel1 );
  thePionMinusInelastic->RegisterMe( theBertini1 );
  thePionMinusInelastic->AddDataSet( thePiXS );
  aProcMan = G4PionMinus::PionMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( thePionMinusInelastic );

  thePionPlusInelastic = new G4PionPlusInelasticProcess();
  thePionPlusInelastic->RegisterMe( theModel1 );
  thePionPlusInelastic->RegisterMe( theBertini1 );
  thePionPlusInelastic->AddDataSet( thePiXS );
  aProcMan = G4PionPlus::PionPlus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( thePionPlusInelastic );

  theKaonMinusInelastic = new G4KaonMinusInelasticProcess();
  theKaonMinusInelastic->RegisterMe( theModel1 );
  theKaonMinusInelastic->RegisterMe( theBertini1 );
  theKaonMinusInelastic->AddDataSet( theKaonXS );
  aProcMan = G4KaonMinus::KaonMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theKaonMinusInelastic );

  theKaonPlusInelastic = new G4KaonPlusInelasticProcess();
  theKaonPlusInelastic->RegisterMe( theModel1 );
  theKaonPlusInelastic->RegisterMe( theBertini1 );
  theKaonPlusInelastic->AddDataSet( theKaonXS );
  aProcMan = G4KaonPlus::KaonPlus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theKaonPlusInelastic );

  theKaonZeroLInelastic = new G4KaonZeroLInelasticProcess();
  theKaonZeroLInelastic->RegisterMe( theModel1 );
  theKaonZeroLInelastic->RegisterMe( theBertini1 );
  theKaonZeroLInelastic->AddDataSet( theKaonXS );
  aProcMan = G4KaonZeroLong::KaonZeroLong()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theKaonZeroLInelastic );

  theKaonZeroSInelastic = new G4KaonZeroSInelasticProcess();
  theKaonZeroSInelastic->RegisterMe( theModel1 );
  theKaonZeroSInelastic->RegisterMe( theBertini1 );
  theKaonZeroSInelastic->AddDataSet( theKaonXS );
  aProcMan = G4KaonZeroShort::KaonZeroShort()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theKaonZeroSInelastic );

  theLambdaInelastic = new G4LambdaInelasticProcess();
  theLambdaInelastic->RegisterMe( theModel2 );
  theLambdaInelastic->RegisterMe( theBertini2 );
  theLambdaInelastic->AddDataSet( theChipsHyperonInelasticXS );
  aProcMan = G4Lambda::Lambda()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theLambdaInelastic );

  theAntiLambdaInelastic = new G4AntiLambdaInelasticProcess();
  theAntiLambdaInelastic->RegisterMe( theModel3 );
  theAntiLambdaInelastic->AddDataSet( theChipsHyperonInelasticXS );
  aProcMan = G4AntiLambda::AntiLambda()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiLambdaInelastic );

  theSigmaMinusInelastic = new G4SigmaMinusInelasticProcess();
  theSigmaMinusInelastic->RegisterMe( theModel2 );
  theSigmaMinusInelastic->RegisterMe( theBertini2 );
  theSigmaMinusInelastic->AddDataSet( theChipsHyperonInelasticXS );
  aProcMan = G4SigmaMinus::SigmaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theSigmaMinusInelastic );

  theAntiSigmaMinusInelastic = new G4AntiSigmaMinusInelasticProcess();
  theAntiSigmaMinusInelastic->RegisterMe( theModel3 );
  theAntiSigmaMinusInelastic->AddDataSet( theChipsHyperonInelasticXS );
  aProcMan = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiSigmaMinusInelastic );

  theSigmaPlusInelastic = new G4SigmaPlusInelasticProcess();
  theSigmaPlusInelastic->RegisterMe( theModel2 );
  theSigmaPlusInelastic->RegisterMe( theBertini2 );
  theSigmaPlusInelastic->AddDataSet( theChipsHyperonInelasticXS );
  aProcMan = G4SigmaPlus::SigmaPlus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theSigmaPlusInelastic );

  theAntiSigmaPlusInelastic = new G4AntiSigmaPlusInelasticProcess();
  theAntiSigmaPlusInelastic->RegisterMe( theModel3 );
  theAntiSigmaPlusInelastic->AddDataSet( theChipsHyperonInelasticXS );
  aProcMan = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiSigmaPlusInelastic );

  theXiMinusInelastic = new G4XiMinusInelasticProcess();
  theXiMinusInelastic->RegisterMe( theModel2 );
  theXiMinusInelastic->RegisterMe( theBertini2 );
  theXiMinusInelastic->AddDataSet( theChipsHyperonInelasticXS );
  aProcMan = G4XiMinus::XiMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theXiMinusInelastic );

  theAntiXiMinusInelastic = new G4AntiXiMinusInelasticProcess();
  theAntiXiMinusInelastic->RegisterMe( theModel3 );
  theAntiXiMinusInelastic->AddDataSet( theChipsHyperonInelasticXS );
  aProcMan = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiXiMinusInelastic );

  theXiZeroInelastic = new G4XiZeroInelasticProcess();
  theXiZeroInelastic->RegisterMe( theModel2 );
  theXiZeroInelastic->RegisterMe( theBertini2 );
  theXiZeroInelastic->AddDataSet( theChipsHyperonInelasticXS );
  aProcMan = G4XiZero::XiZero()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theXiZeroInelastic );

  theAntiXiZeroInelastic = new G4AntiXiZeroInelasticProcess();
  theAntiXiZeroInelastic->RegisterMe( theModel3 );
  theAntiXiZeroInelastic->AddDataSet( theChipsHyperonInelasticXS );
  aProcMan = G4AntiXiZero::AntiXiZero()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiXiZeroInelastic );

  theOmegaMinusInelastic = new G4OmegaMinusInelasticProcess();
  theOmegaMinusInelastic->RegisterMe( theModel2 );
  theOmegaMinusInelastic->RegisterMe( theBertini2 );
  theOmegaMinusInelastic->AddDataSet( theChipsHyperonInelasticXS );
  aProcMan = G4OmegaMinus::OmegaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theOmegaMinusInelastic );

  theAntiOmegaMinusInelastic = new G4AntiOmegaMinusInelasticProcess();
  theAntiOmegaMinusInelastic->RegisterMe( theModel3 );
  theAntiOmegaMinusInelastic->AddDataSet( theChipsHyperonInelasticXS );
  aProcMan = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiOmegaMinusInelastic );

  theAntiProtonInelastic = new G4AntiProtonInelasticProcess();
  theAntiProtonInelastic->RegisterMe( theModel3 );
  theAntiProtonInelastic->AddDataSet( theAntiNucleonXS );
  aProcMan = G4AntiProton::AntiProton()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiProtonInelastic );

  theAntiNeutronInelastic = new G4AntiNeutronInelasticProcess();
  theAntiNeutronInelastic->RegisterMe( theModel3 );
  theAntiNeutronInelastic->AddDataSet( theAntiNucleonXS );
  aProcMan = G4AntiNeutron::AntiNeutron()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiNeutronInelastic );

  theAntiDeuteronInelastic = new G4AntiDeuteronInelasticProcess();
  theAntiDeuteronInelastic->RegisterMe( theModel3 );
  theAntiDeuteronInelastic->AddDataSet( theAntiNucleonXS );
  aProcMan = G4AntiDeuteron::AntiDeuteron()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiDeuteronInelastic );

  theAntiTritonInelastic = new G4AntiTritonInelasticProcess();
  theAntiTritonInelastic->RegisterMe( theModel3 );
  theAntiTritonInelastic->AddDataSet( theAntiNucleonXS );
  aProcMan = G4AntiTriton::AntiTriton()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiTritonInelastic );

  theAntiHe3Inelastic = new G4AntiHe3InelasticProcess();
  theAntiHe3Inelastic->RegisterMe( theModel3 );
  theAntiHe3Inelastic->AddDataSet( theAntiNucleonXS );
  aProcMan = G4AntiHe3::AntiHe3()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiHe3Inelastic );

  theAntiAlphaInelastic = new G4AntiAlphaInelasticProcess();
  theAntiAlphaInelastic->RegisterMe( theModel3 );
  theAntiAlphaInelastic->AddDataSet( theAntiNucleonXS );
  aProcMan = G4AntiAlpha::AntiAlpha()->GetProcessManager();
  aProcMan->AddDiscreteProcess( theAntiAlphaInelastic );
}

