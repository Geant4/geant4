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
#include "G4CrossSectionInelastic.hh"
#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4HadronicParameters.hh"
  
  
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
  
  HyperonFTFP->SetMinEnergy( G4HadronicParameters::Instance()->GetMinEnergyTransitionFTF_Cascade() );
  HyperonFTFP->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );

  theStringModel = new G4FTFModel;
  theStringDecay = new G4ExcitedStringDecay(theLund = new G4LundStringFragmentation);
  theStringModel->SetFragmentationModel(theStringDecay);

  theCascade = new G4GeneratorPrecompoundInterface;

  HyperonFTFP->SetTransport(theCascade);
  HyperonFTFP->SetHighEnergyGenerator(theStringModel);
  
  theBertini = new G4CascadeInterface;
  theBertini->SetMinEnergy( 0.0 );
  theBertini->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergyTransitionFTF_Cascade() );

  // AntiHyperons: Use FTFP for full energy range, starting at 0.  

  AntiHyperonFTFP = new G4TheoFSGenerator("FTFP");
  AntiHyperonFTFP->SetMinEnergy( 0.0 );
  AntiHyperonFTFP->SetMaxEnergy( G4HadronicParameters::Instance()->GetMaxEnergy() );
  AntiHyperonFTFP->SetTransport(theCascade);
  AntiHyperonFTFP->SetHighEnergyGenerator(theStringModel);

  // use Glauber-Gribov cross sections
  theInelasticCrossSection = new G4CrossSectionInelastic( new G4ComponentGGHadronNucleusXsc );
}

G4HyperonFTFPBuilder::~G4HyperonFTFPBuilder()
{
  delete theStringDecay;
  delete theLund;
}

void G4HyperonFTFPBuilder::Build()
{
  G4ProcessManager * aProcMan = 0;
  wasActivated = true;

  // Lambda
  theLambdaInelastic = new G4LambdaInelasticProcess();
  theLambdaInelastic->RegisterMe(theBertini);
  theLambdaInelastic->RegisterMe(HyperonFTFP);
  theLambdaInelastic->AddDataSet(theInelasticCrossSection);
  aProcMan = G4Lambda::Lambda()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theLambdaInelastic);
  
  // AntiLambda
  theAntiLambdaInelastic = new G4AntiLambdaInelasticProcess();
  theAntiLambdaInelastic->RegisterMe(AntiHyperonFTFP);
  theAntiLambdaInelastic->AddDataSet(theInelasticCrossSection);
  
  aProcMan = G4AntiLambda::AntiLambda()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theAntiLambdaInelastic);
    
  // SigmaMinus
  theSigmaMinusInelastic = new G4SigmaMinusInelasticProcess();
  theSigmaMinusInelastic->RegisterMe(theBertini);
  theSigmaMinusInelastic->RegisterMe(HyperonFTFP);
  theSigmaMinusInelastic->AddDataSet(theInelasticCrossSection);

  aProcMan = G4SigmaMinus::SigmaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theSigmaMinusInelastic);

  // anti-SigmaMinus
  theAntiSigmaMinusInelastic = new G4AntiSigmaMinusInelasticProcess();
  theAntiSigmaMinusInelastic->RegisterMe(AntiHyperonFTFP);
  theAntiSigmaMinusInelastic->AddDataSet(theInelasticCrossSection);

  aProcMan = G4AntiSigmaMinus::AntiSigmaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theAntiSigmaMinusInelastic);

  // SigmaPlus
  theSigmaPlusInelastic = new G4SigmaPlusInelasticProcess();
  theSigmaPlusInelastic->RegisterMe(theBertini);
  theSigmaPlusInelastic->RegisterMe(HyperonFTFP);
  theSigmaPlusInelastic->AddDataSet(theInelasticCrossSection);

  aProcMan = G4SigmaPlus::SigmaPlus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theSigmaPlusInelastic);

  // anti-SigmaPlus
  theAntiSigmaPlusInelastic = new G4AntiSigmaPlusInelasticProcess();
  theAntiSigmaPlusInelastic->RegisterMe(AntiHyperonFTFP);
  theAntiSigmaPlusInelastic->AddDataSet(theInelasticCrossSection);

  aProcMan = G4AntiSigmaPlus::AntiSigmaPlus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theAntiSigmaPlusInelastic);

  // XiMinus
  theXiMinusInelastic = new G4XiMinusInelasticProcess();
  theXiMinusInelastic->RegisterMe(theBertini);
  theXiMinusInelastic->RegisterMe(HyperonFTFP);
  theXiMinusInelastic->AddDataSet(theInelasticCrossSection);

  aProcMan = G4XiMinus::XiMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theXiMinusInelastic);

  // anti-XiMinus
  theAntiXiMinusInelastic = new G4AntiXiMinusInelasticProcess();
  theAntiXiMinusInelastic->RegisterMe(AntiHyperonFTFP);
  theAntiXiMinusInelastic->AddDataSet(theInelasticCrossSection);

  aProcMan = G4AntiXiMinus::AntiXiMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theAntiXiMinusInelastic);

  // XiZero
  theXiZeroInelastic = new G4XiZeroInelasticProcess();
  theXiZeroInelastic->RegisterMe(theBertini);
  theXiZeroInelastic->RegisterMe(HyperonFTFP);
  theXiZeroInelastic->AddDataSet(theInelasticCrossSection);

  aProcMan = G4XiZero::XiZero()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theXiZeroInelastic);

  // anti-XiZero
  theAntiXiZeroInelastic = new G4AntiXiZeroInelasticProcess();
  theAntiXiZeroInelastic->RegisterMe(AntiHyperonFTFP);
  theAntiXiZeroInelastic->AddDataSet(theInelasticCrossSection);

  aProcMan = G4AntiXiZero::AntiXiZero()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theAntiXiZeroInelastic);

  // OmegaMinus
  theOmegaMinusInelastic = new G4OmegaMinusInelasticProcess();
  theOmegaMinusInelastic->RegisterMe(theBertini);
  theOmegaMinusInelastic->RegisterMe(HyperonFTFP);
  theOmegaMinusInelastic->AddDataSet(theInelasticCrossSection);

  aProcMan = G4OmegaMinus::OmegaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theOmegaMinusInelastic);

  // anti-OmegaMinus
  theAntiOmegaMinusInelastic = new G4AntiOmegaMinusInelasticProcess();
  theAntiOmegaMinusInelastic->RegisterMe(AntiHyperonFTFP);
  theAntiOmegaMinusInelastic->AddDataSet(theInelasticCrossSection);

  aProcMan = G4AntiOmegaMinus::AntiOmegaMinus()->GetProcessManager();
  aProcMan->AddDiscreteProcess(theAntiOmegaMinusInelastic);
}

