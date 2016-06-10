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
//---------------------------------------------------------------------------
//
// ClassName:   G4HyperonFTFPBuilder
//
// Author: 2012 G.Folger
//    Implementation started from G4HyperonLHEPBuilder.  
//
// Modified:
//----------------------------------------------------------------------------
//
#include "G4HyperonFTFPBuilder.hh"

#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4CrossSectionDataSetRegistry.hh"
  
  
G4HyperonFTFPBuilder::G4HyperonFTFPBuilder(): 
 theLambdaInelastic(0),
 theAntiLambdaInelastic(0),
 theSigmaMinusInelastic(0),
 theAntiSigmaMinusInelastic(0),
 theSigmaPlusInelastic(0), 
 theAntiSigmaPlusInelastic(0), 
 theXiZeroInelastic(0), 
 theAntiXiZeroInelastic(0),
 theXiMinusInelastic(0), 
 theAntiXiMinusInelastic(0),
 theOmegaMinusInelastic(0), 
 theAntiOmegaMinusInelastic(0), 
 wasActivated(false)
{

// Hyperon : Bertini at low energies, then FTFP

  HyperonFTFP = new G4TheoFSGenerator("FTFP");
  
  HyperonFTFP->SetMinEnergy( 2.*GeV );
  HyperonFTFP->SetMaxEnergy( 100.*TeV );

  theStringModel = new G4FTFModel;
  theStringDecay = new G4ExcitedStringDecay(theLund = new G4LundStringFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);

  theCascade = new G4GeneratorPrecompoundInterface;

  HyperonFTFP->SetTransport(theCascade);
  HyperonFTFP->SetHighEnergyGenerator(theStringModel);
  
  theBertini = new G4CascadeInterface;
  theBertini->SetMinEnergy( 0.*GeV );
  theBertini->SetMaxEnergy( 6.*GeV );

// AntiHyperons: Use FTFP for full energy range, starting at 0.  

  AntiHyperonFTFP = new G4TheoFSGenerator("FTFP");
  AntiHyperonFTFP->SetMinEnergy( 0.*GeV );
  AntiHyperonFTFP->SetMaxEnergy( 100.*TeV );
  AntiHyperonFTFP->SetTransport(theCascade);
  AntiHyperonFTFP->SetHighEnergyGenerator(theStringModel);

// use CHIPS cross sections
  theCHIPSInelastic = G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsHyperonInelasticXS::Default_Name());
}


G4HyperonFTFPBuilder::~G4HyperonFTFPBuilder()
{
  //delete HyperonFTFP;
  delete theStringModel;
  delete theStringDecay;
  delete theLund;
  //delete AntiHyperonFTFP;
  /*  
  if (wasActivated) {
     delete theLambdaInelastic;
     delete theAntiLambdaInelastic;
     delete theSigmaMinusInelastic;
     delete theAntiSigmaMinusInelastic;
     delete theSigmaPlusInelastic;
     delete theAntiSigmaPlusInelastic;
     delete theXiMinusInelastic;
     delete theAntiXiMinusInelastic;
     delete theXiZeroInelastic;
     delete theAntiXiZeroInelastic;
     delete theOmegaMinusInelastic;
     delete theAntiOmegaMinusInelastic;
  } 
  */  
}

void G4HyperonFTFPBuilder::Build()
{
  G4ProcessManager * aProcMan = 0;
  wasActivated = true;

  // Lambda
  theLambdaInelastic = new G4LambdaInelasticProcess();
  theLambdaInelastic->RegisterMe(theBertini);
  theLambdaInelastic->RegisterMe(HyperonFTFP);
  theLambdaInelastic->AddDataSet(theCHIPSInelastic);
  aProcMan = G4Lambda::Lambda()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theLambdaInelastic);
  
  // AntiLambda
  theAntiLambdaInelastic = new G4AntiLambdaInelasticProcess();
  theAntiLambdaInelastic->RegisterMe(AntiHyperonFTFP);
  theAntiLambdaInelastic->AddDataSet(theCHIPSInelastic);
  
  aProcMan = G4AntiLambda::AntiLambda()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theAntiLambdaInelastic);
    
  // SigmaMinus
  theSigmaMinusInelastic = new G4SigmaMinusInelasticProcess();
  theSigmaMinusInelastic->RegisterMe(theBertini);
  theSigmaMinusInelastic->RegisterMe(HyperonFTFP);
  theSigmaMinusInelastic->AddDataSet(theCHIPSInelastic);

  aProcMan = G4SigmaMinus::SigmaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theSigmaMinusInelastic);

  // anti-SigmaMinus
  theAntiSigmaMinusInelastic = new G4AntiSigmaMinusInelasticProcess();
  theAntiSigmaMinusInelastic->RegisterMe(AntiHyperonFTFP);
  theAntiSigmaMinusInelastic->AddDataSet(theCHIPSInelastic);

  aProcMan = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theAntiSigmaMinusInelastic);

  // SigmaPlus
  theSigmaPlusInelastic = new G4SigmaPlusInelasticProcess();
  theSigmaPlusInelastic->RegisterMe(theBertini);
  theSigmaPlusInelastic->RegisterMe(HyperonFTFP);
  theSigmaPlusInelastic->AddDataSet(theCHIPSInelastic);

  aProcMan = G4SigmaPlus::SigmaPlus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theSigmaPlusInelastic);

  // anti-SigmaPlus
  theAntiSigmaPlusInelastic = new G4AntiSigmaPlusInelasticProcess();
  theAntiSigmaPlusInelastic->RegisterMe(AntiHyperonFTFP);
  theAntiSigmaPlusInelastic->AddDataSet(theCHIPSInelastic);

  aProcMan = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theAntiSigmaPlusInelastic);

  // XiMinus
  theXiMinusInelastic = new G4XiMinusInelasticProcess();
  theXiMinusInelastic->RegisterMe(theBertini);
  theXiMinusInelastic->RegisterMe(HyperonFTFP);
  theXiMinusInelastic->AddDataSet(theCHIPSInelastic);

  aProcMan = G4XiMinus::XiMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theXiMinusInelastic);

  // anti-XiMinus
  theAntiXiMinusInelastic = new G4AntiXiMinusInelasticProcess();
  theAntiXiMinusInelastic->RegisterMe(AntiHyperonFTFP);
  theAntiXiMinusInelastic->AddDataSet(theCHIPSInelastic);

  aProcMan = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theAntiXiMinusInelastic);

  // XiZero
  theXiZeroInelastic = new G4XiZeroInelasticProcess();
  theXiZeroInelastic->RegisterMe(theBertini);
  theXiZeroInelastic->RegisterMe(HyperonFTFP);
  theXiZeroInelastic->AddDataSet(theCHIPSInelastic);

  aProcMan = G4XiZero::XiZero()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theXiZeroInelastic);

  // anti-XiZero
  theAntiXiZeroInelastic = new G4AntiXiZeroInelasticProcess();
  theAntiXiZeroInelastic->RegisterMe(AntiHyperonFTFP);
  theAntiXiZeroInelastic->AddDataSet(theCHIPSInelastic);

  aProcMan = G4AntiXiZero::AntiXiZero()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theAntiXiZeroInelastic);

  // OmegaMinus
  theOmegaMinusInelastic = new G4OmegaMinusInelasticProcess();
  theOmegaMinusInelastic->RegisterMe(theBertini);
  theOmegaMinusInelastic->RegisterMe(HyperonFTFP);
  theOmegaMinusInelastic->AddDataSet(theCHIPSInelastic);

  aProcMan = G4OmegaMinus::OmegaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theOmegaMinusInelastic);

  // anti-OmegaMinus
  theAntiOmegaMinusInelastic = new G4AntiOmegaMinusInelasticProcess();
  theAntiOmegaMinusInelastic->RegisterMe(AntiHyperonFTFP);
  theAntiOmegaMinusInelastic->AddDataSet(theCHIPSInelastic);

  aProcMan = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theAntiOmegaMinusInelastic);
}


